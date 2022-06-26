#pragma once
#include <Engine/Scene/Entity.hpp>
#include <Eigen/src/Core/AssignEvaluator.h>
#include <Eigen/src/Core/GenericPacketMath.h>
#include <Core/CoreMacros.hpp>
#include <Core/Types.hpp>
#include <Core/Utils/Observable.hpp>
#include <mutex>
#include <string>
#include <vector>

namespace Ra {
namespace Engine {
namespace Scene {
class Entity;

inline const std::string& Entity::getName() const {
    return m_name;
}

inline void Entity::rename( const std::string& name ) {
    m_name = name;
}

inline void Entity::setTransform( const Core::Transform& transform ) {
    m_transformChanged        = true;
    m_doubleBufferedTransform = transform;
}

inline void Entity::setTransform( const Core::Matrix4& transform ) {
    setTransform( Core::Transform( transform ) );
}

inline const Core::Transform& Entity::getTransform() const {
    // Radium-V2 : why a mutex on read ? there is no lock on write on this!
    std::lock_guard<std::mutex> lock( m_transformMutex );
    return m_transform;
}

inline const Core::Matrix4& Entity::getTransformAsMatrix() const {
    // Radium-V2 : why a mutex on read ? there is no lock on write on this!
    std::lock_guard<std::mutex> lock( m_transformMutex );
    return m_transform.matrix();
}

inline uint Entity::getNumComponents() const {
    return uint( m_components.size() );
}

inline Core::Utils::Observable<const Entity*>& Entity::transformationObservers() const {
    return m_transformationObservers;
}

} // namespace Scene
} // namespace Engine
} // namespace Ra
