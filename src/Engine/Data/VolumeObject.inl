#pragma once
#include <Engine/Data/VolumeObject.hpp>
#include <Core/CoreMacros.hpp>
#include <Core/Geometry/AbstractGeometry.hpp>
#include <memory>

namespace Ra {
namespace Core {
namespace Geometry {
class AbstractVolume;
}  // namespace Geometry
}  // namespace Core

namespace Engine {
namespace Data {
const Core::Geometry::AbstractGeometry& VolumeObject::getAbstractGeometry() const {
    CORE_ASSERT( m_volume, "Volume is not initialized" );
    return *(Core::Geometry::AbstractGeometry*)( m_volume.get() );
}

Core::Geometry::AbstractGeometry& VolumeObject::getAbstractGeometry() {
    CORE_ASSERT( m_volume, "Volume is not initialized" );
    return *(Core::Geometry::AbstractGeometry*)( m_volume.get() );
}

/// Returns the underlying AbstractVolume
const Core::Geometry::AbstractVolume& VolumeObject::getVolume() const {
    CORE_ASSERT( m_volume, "Volume is not initialized" );
    return *m_volume.get();
}

Core::Geometry::AbstractVolume& VolumeObject::getVolume() {
    CORE_ASSERT( m_volume, "Volume is not initialized" );
    return *m_volume.get();
}

} // namespace Data
} // namespace Engine
} // namespace Ra
