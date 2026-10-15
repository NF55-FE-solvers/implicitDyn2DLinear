#include <iostream>
#include <fstream>
#include "Eigen/Eigen"
#include "Eigen/PardisoSupport"
#include "Mesh.h"
#include "import.h"
#include "settings.h"
#include "TimerHighRes.h"
#include "estimateMemoryUsage.h"
#include "material.h"
#include "Bc.h"
#include "rearrangeDof.h"
#include "elMatrices.h"
#include "export.h"
#include "Ic.h"
#include "clearResults.h"

int main(){
    //Timer timer {};
    //timer.start();
    Mesh mesh {};
    Bc bc {}; // boundary conditions class
    using ndofVectf = Eigen::Vector<Settings::precision, Eigen::Dynamic>;
    using ndofVecti = Eigen::Vector<int, Eigen::Dynamic>;
    using SizT = std::size_t;
    ndofVectf dispOld;
    ndofVectf force;
    
    // To handle exceptions in any function/method
    try{
        //-----Importing mesh and boundary conditions
        importMesh(mesh, bc);
        std::cout << "Problem size: " << mesh.nDof << " d.o.f.\n";
        double estimatedMemoryUsageMB {estimateMemUsageMB(mesh.nDof
            , mesh.nElements)};
        std::cout << "Esitmated Memory usage (before solving phase): " 
            << estimatedMemoryUsageMB/Settings::multiplierMemory4Solving 
            << " MB\nEstimated MAX Memory usage: " << estimatedMemoryUsageMB 
            << " MB\n\n";
        if (estimatedMemoryUsageMB > Settings::memoryLimitMB)
            throw "Estimated memory usage > memory usage limit set in settings";

        //-----(global mass matrix + global stiffness matrix) computation
        Eigen::SparseMatrix<Settings::precision> kG (mesh.nDof, mesh.nDof);
        Eigen::SparseMatrix<Settings::precision> mG (mesh.nDof, mesh.nDof);
        {
            std::vector<Eigen::Triplet<Settings::precision>> tripletListK;
            std::vector<Eigen::Triplet<Settings::precision>> tripletListM;
            tripletListK.reserve(mesh.nElements*Settings::nDofKel);
            tripletListM.reserve(mesh.nElements*Settings::nDofKel);
            kmMat::KMelEigen kEl;
            kmMat::KMelEigen mEl;
            kmMat::NodesCoordEigen nodesCoords;
            kmMat::NEigen N;
            kmMat::dNEigen dN;
            kmMat::invJEigen invJ;
            kmMat::BEigen bEl;
            Eigen::Matrix<int,Settings::dofNode,Settings::nNodesEl> 
                dofEliMat;
            int kmElContiguousIndx {};
            const float numEl1perc {mesh.nElements/20.0f};
            float numEl2Message {numEl1perc};
            for (SizT e=0uz; e<static_cast<SizT>(mesh.nElements); ++e){
                nodesCoords = mesh.nodes(mesh.elements.row(e)
                    ,Eigen::placeholders::all);
                kmMat::computeElKM(nodesCoords,N,dN,invJ,bEl,kEl,mEl);
                dofEliMat.row(0) = mesh.elements.row(e)*2;
                dofEliMat.row(1) = mesh.elements.row(e)*2 
                    + Eigen::RowVector<int,Settings::nNodesEl>::Ones();
                const auto dofEli = dofEliMat.reshaped(); // view
                // as Eigen Matrices are column-major, stored contiguously 
                // column by column, it's more efficient to read them also 
                // in this way
                kmElContiguousIndx = 0;
                for (const int j : dofEli){
                    for (const int i : dofEli){
                        tripletListM.emplace_back(i,j
                            ,mEl(kmElContiguousIndx));
                        tripletListK.emplace_back(i,j
                            ,kEl(kmElContiguousIndx));
                        ++kmElContiguousIndx;
                    }
                }
                if (e > numEl2Message){
                    std::cout << "Sparse global mass & stiffness matrix " 
                        << "computation... " << numEl2Message/numEl1perc
                        *2.0f << "%\n";
                    numEl2Message += numEl1perc;
                }
            }
            kG.setFromTriplets(tripletListK.begin(), tripletListK.end());
            mG.setFromTriplets(tripletListM.begin(), tripletListM.end());
            std::cout << "\nGlobal Mass [M] and Stiffness [K] matrices computed"
                << " successfully!\n\n";
        } // Triplet list and all other helper matrices deallocated

        //-----Solving system of eqs. (splitting method) 
        //     & time integration (Newmark-beta method)

        dispOld = ndofVectf::Constant(mesh.nDof,Settings::NaN);
        force = ndofVectf::Zero(mesh.nDof);

        bc.setBcInit(dispOld, force);

        // Re-arranging vectors (splitting method) + saving mapping vector
        ndofVecti mapDof;
        mapDof = ndofVecti::LinSpaced(mesh.nDof, 0, mesh.nDof-1);
        rearrange(dispOld, force, mapDof);
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> 
            mapDofMat {mapDof};
        kG = kG.twistedBy(mapDofMat.transpose()); // Re-arranging
        mG = mG.twistedBy(mapDofMat.transpose()); // Re-arranging

        int lastNanIndx {mesh.nDof-1};
        for (; lastNanIndx>0; --lastNanIndx){
            if (std::isnan(dispOld(lastNanIndx))) break;
        }
        int nUnknownDisp {lastNanIndx+1};
        int nKnownDisp {mesh.nDof - nUnknownDisp};

        ndofVectf dispNew, velOld, velNew, accOld, accNew;

        Settings::solverType solver;
        
        Ic::setIC(dispOld, velOld, accOld, force, mG, kG, solver, nUnknownDisp);
        dispNew = dispOld;
        velNew  = velOld;
        accNew  = accOld;

        Settings::precision time {Settings::t0 + Settings::dt};
        Settings::precision tSaveVTK {Settings::t0 + Settings::dtSaveVTK};

        clearPriorResults();

        while (time < (Settings::tFinal + Settings::dt)){
            

            time += Settings::dt;
        }

        // 6-Change the Bc.cpp(+h) so to be able to update the BC/loads in time
        // 7-Solve for the current time instant
        // 8-Write the results in a .vtk file
        // 9-Update the variables with std::swap(old,new) + advance in time


        solver.compute(kG.block(0,0,lastNanIndx+1,lastNanIndx+1));
        if(solver.info()!=Eigen::Success) {
            throw "Linear solver failed decomposition(1).";
        }
        Eigen::Vector<Settings::precision,Eigen::Dynamic> rhsDispSol = 
            force.head(nUnknownDisp);
        rhsDispSol.noalias() -= kG.topRightCorner(nUnknownDisp,nKnownDisp)
            *dispOld.tail(nKnownDisp);
        dispOld.head(nUnknownDisp) = solver.solve(rhsDispSol);
        if(solver.info()!=Eigen::Success) {
            throw "Linear solver failed solving(1).";
        }
        std::cout << "Unknown nodal displacements computed successfully!"
            << "\n\n";
        dispOld = mapDofMat * dispOld;

        //-----Exporting the results
        if (Settings::exportBinaryVTK) Export::exportResultsBinary(mesh, dispOld);
        else Export::exportResultsASCII(mesh, dispOld);
        
        //timer.stop();
        //std::cout << "Elapsed: " << std::setprecision(10) 
        //    << timer.getDurationMicro() << " micro sec.\n";
    }
    catch(const char* exceptionMsg){
        std::cerr << "Error: " << exceptionMsg << '\n';
    }
    catch(...){
        std::cerr << "Unrecognized error occured!\n";
    }

    return 0;
}
