#pragma once
#include <Engine/Data/RawShaderMaterial.hpp>
#include <Engine/Data/RenderParameters.hpp>
#include <memory>

namespace Ra {
namespace Core {
namespace Asset {
class MaterialData;
class RawShaderMaterialData;
}  // namespace Asset
}  // namespace Core

namespace Engine {
namespace Data {
class Material;

inline Material*
RawShaderMaterialConverter::operator()( const Ra::Core::Asset::MaterialData* toconvert ) {
    auto mat = static_cast<const Core::Asset::RawShaderMaterialData*>( toconvert );
    return new RawShaderMaterial( mat->getName(), mat->m_shaders, mat->m_paramProvider );
}

inline Data::RenderParameters& RawShaderMaterial::getParameters() {
    return m_paramProvider->getParameters();
}

inline const Data::RenderParameters& RawShaderMaterial::getParameters() const {
    return m_paramProvider->getParameters();
}

} // namespace Data
} // namespace Engine
} // namespace Ra
