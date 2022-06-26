#pragma once

#include <Core/Asset/Camera.hpp>
#include <Engine/RaEngine.hpp>
#include <Engine/Scene/CameraComponent.hpp>
#include <Engine/Scene/CameraManager.hpp>
#include <Engine/Scene/CameraStorage.hpp>
#include <memory>
#include <stddef.h>
#include <vector>

namespace Ra {
namespace Engine {
namespace Scene {
class CameraComponent;

/**
 * Associated class.
 */
class RA_ENGINE_API DefaultCameraStorage : public CameraStorage
{
  public:
    DefaultCameraStorage();
    void add( CameraComponent* cam ) override;
    void remove( CameraComponent* cam ) override;
    size_t size() const override;
    void clear() override;
    CameraComponent* operator[]( unsigned int n ) override;

  private:
    /** Vectors of Camera references. */
    // std::multmmap<Ra::Core::Asset::Camera::ProjType, CameraComponent*> m_data;
    std::vector<CameraComponent*> m_data;
};

/**
 * @brief DefaultCameraManager. A simple Camera Manager with a list of Cameras.
 */
class RA_ENGINE_API DefaultCameraManager : public CameraManager
{
  public:
    DefaultCameraManager();

    /// Return the \p cam-th camera.
    const CameraComponent* getCamera( size_t cam ) const override;

    /// Add \p cam for management.
    void addCamera( CameraComponent* cam ) override;
};

} // namespace Scene
} // namespace Engine
} // namespace Ra
