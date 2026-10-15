#ifndef SETTINGS_H
#define SETTINGS_H

#include <limits>
#include <Eigen/Core>
#include <filesystem>

namespace Settings{
    // main aliases
    using precision = double; // for computations (float or double), leave
    // double when using Intel MKL Pardiso solvers
    // Suggested for square positive definite system:
    // For problems up to 1M dof -> PardisoLLT (direct)
    // For large problems (over 1-2M dof) -> ConjugateGradient (iterative)
    using solverType = Eigen::PardisoLLT<Eigen::SparseMatrix<precision>>;
    using fPath = std::filesystem::path;

    // for time integration
    constexpr const precision t0 {0.0};
    constexpr const precision tFinal {1.0};
    constexpr const precision dt {1e-2};
    constexpr const int timeInst2save {50};
    constexpr const precision dtSaveVTK {tFinal/timeInst2save};
    constexpr const precision gamma {0.5}; // const. of Newmark-beta method
    constexpr const precision beta {0.25}; // const. of Newmark-beta method
    
    // for results and other things
    constexpr const bool exportBinaryVTK{true}; // binary is much faster
    constexpr const bool fullPrecisionASCIIWriting {false}; // suggested false
    constexpr const double memoryLimitMB {3000.0}; // avoiding dynamic 
    // allocation errors due to lack of heap memory availability
    fPath expectedRootProj {"implicitDyn2DLinear"};
    fPath resultsDirName {"./VTK RESULTS"};


    // Useful constants:
    constexpr const double multiplierMemory4Solving {2.5}; // the multiplier of 
    // memory estimated before solving, to account for the one of solving phase
    constexpr const precision NaN {std::numeric_limits<precision>::quiet_NaN()};
    constexpr const precision sqrt3 {1.73205080756887729352};
    constexpr const std::size_t nNodesEl {4uz};
    constexpr const std::size_t dofNode {2uz};
    constexpr const std::size_t nDofEl {nNodesEl * dofNode};
    constexpr const std::size_t nDofKel {nDofEl * nDofEl};
}

#endif