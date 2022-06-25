#include <Core/Geometry/HeatDiffusion.hpp>
#include <Eigen/src/Core/Assign.h>
#include <Eigen/src/Core/AssignEvaluator.h>
#include <Eigen/src/Core/CwiseBinaryOp.h>
#include <Eigen/src/Core/CwiseNullaryOp.h>
#include <Eigen/src/Core/DenseCoeffsBase.h>
#include <Eigen/src/Core/DiagonalMatrix.h>
#include <Eigen/src/Core/GenericPacketMath.h>
#include <Eigen/src/Core/MathFunctions.h>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/PermutationMatrix.h>
#include <Eigen/src/Core/Redux.h>
#include <Eigen/src/Core/Solve.h>
#include <Eigen/src/Core/arch/SSE/PacketMath.h>
#include <Eigen/src/Core/functors/BinaryFunctors.h>
#include <Eigen/src/Core/util/IntegralConstant.h>
#include <Eigen/src/Core/util/Memory.h>
#include <Eigen/src/Core/util/XprHelper.h>
#include <Eigen/src/OrderingMethods/Amd.h>
#include <Eigen/src/SparseCholesky/SimplicialCholesky.h>
#include <Eigen/src/SparseCore/SparseAssign.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <Eigen/src/SparseCore/SparseMatrixBase.h>
#include <Eigen/src/SparseCore/SparseSelfAdjointView.h>
#include <Eigen/src/SparseCore/TriangularSolver.h>
#include <algorithm>
#include <new>
#include <string.h>
#include <utility>

#include "Core/Containers/VectorArray.hpp"
#include "Core/Geometry/Area.hpp"
#include "Core/Geometry/Laplacian.hpp"

namespace Ra {
namespace Core {
namespace Geometry {

Time t( const Scalar& m, const Scalar& h ) {
    return ( m * h * h );
}

void heat( const AreaMatrix& A,
           const Time& t,
           const LaplacianMatrix& L,
           Heat& u,
           const Sparse& delta ) {
    Eigen::SimplicialLLT<Sparse> llt;
    llt.compute( ( A + ( t * L ) ) );
    VectorN b  = delta;
    u.getMap() = llt.solve( b );
}

Heat heat( const AreaMatrix& A, const Time& t, const LaplacianMatrix& L, const Sparse& delta ) {
    Heat u( L.rows() );
    Eigen::SimplicialLLT<Sparse> llt;
    llt.compute( ( A + ( t * L ) ) );
    VectorN b  = delta;
    u.getMap() = llt.solve( b );
    return u;
}

} // namespace Geometry
} // namespace Core
} // namespace Ra
