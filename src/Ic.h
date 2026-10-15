#ifndef IC_H
#define IC_H

#include <iostream>
#include "Eigen/Eigen"
#include "settings.h"

namespace Ic{
    using ndofVectf = Eigen::Vector<Settings::precision, Eigen::Dynamic>;
    using sparseMat = Eigen::SparseMatrix<Settings::precision>;

    #define IMPORT_FROM_STATIC 0

    void setIC(ndofVectf& dispOld, ndofVectf& velOld, ndofVectf& accOld
        , ndofVectf& force, sparseMat& mG, sparseMat& kG
        , Settings::solverType& solver, int nUnkDisp)
    {
        std::cout << "\nEmbedding initial conditions @t=" << Settings::t0 
            << "s ...\n";
        #if IMPORT_FROM_STATIC
            // TBD: after the import of nodal displacements, set velOld and for
            //      accOld compute inv(M)*(F-K*dispOld)
        #else
            dispOld.head(nUnkDisp) = ndofVectf::Zero(nUnkDisp);
            velOld = ndofVectf::Zero(dispOld.size());
            accOld = ndofVectf::Zero(dispOld.size());
        #endif
    }
}

#endif