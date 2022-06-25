#pragma once

#include <Core/Containers/VectorArray.hpp>
#include <Core/Types.hpp>
#include <vector>

#include "Core/CoreMacros.hpp"

namespace Ra {
namespace Core {
template <typename V>
class VectorArray;

namespace Geometry {

// get distance between two sets of vertices.
void vertexDistance( const VectorArray<Vector3>& v0,
                     const VectorArray<Vector3>& v1,
                     std::vector<Scalar>& sqrDist,
                     Scalar& sqrMin,
                     Scalar& sqrMax,
                     Scalar& sqrAvg );

void vertexDistance( const VectorArray<Vector3>& v0,
                     const VectorArray<Vector3>& v1,
                     Scalar& sqrMin,
                     Scalar& sqrMax,
                     Scalar& sqrAvg );

Scalar vertexDistance( const VectorArray<Vector3>& v0, const VectorArray<Vector3>& v1 );

} // namespace Geometry
} // namespace Core
} // namespace Ra
