#pragma once
#include "DirLight.hpp"

#include <Eigen/src/Core/AssignEvaluator.h>
#include <Eigen/src/Core/Dot.h>
#include <Eigen/src/Core/MathFunctions.h>
#include <Eigen/src/Core/Matrix.h>
#include <Core/CoreMacros.hpp>

namespace Ra {
namespace Engine {
namespace Scene {

inline void DirectionalLight::setDirection( const Eigen::Matrix<Scalar, 3, 1>& dir ) {
    m_direction = dir.normalized();
}

inline const Eigen::Matrix<Scalar, 3, 1>& DirectionalLight::getDirection() const {
    return m_direction;
}

} // namespace Scene
} // namespace Engine
} // namespace Ra
