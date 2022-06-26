#include <Core/Containers/VectorArray.hpp>
#include <Core/CoreMacros.hpp>
#include <Core/Geometry/MeshPrimitives.hpp>
#include <Core/Geometry/StandardAttribNames.hpp>
#include <Core/Geometry/TriangleMesh.hpp>
#include <Core/Geometry/Volume.hpp>
#include <Core/Types.hpp>
#include <Core/Utils/Attribs.hpp>
#include <Eigen/src/Core/CommaInitializer.h>
#include <Eigen/src/Core/DenseBase.h>
#include <Eigen/src/Core/DenseCoeffsBase.h>
#include <Eigen/src/Core/Map.h>
#include <Eigen/src/Core/Matrix.h>
#include <Engine/Data/Mesh.hpp>
#include <Engine/Data/VolumeObject.hpp>
#include <Engine/OpenGL.hpp>
#include <glbinding/Boolean8.h>
#include <glbinding/gl/boolean.h>
#include <glbinding/gl/enum.h>
#include <glbinding/gl/functions.h>
#include <glbinding/gl45core/enum.h>
#include <glbinding/gl45core/functions.h>

namespace Ra {
namespace Engine {
namespace Data {
class ShaderProgram;

VolumeObject::VolumeObject( const std::string& name ) :
    Displayable( name ), m_tex( {} ), m_mesh( name + "_internal" ) {}

VolumeObject::~VolumeObject() {}

void VolumeObject::loadGeometry( Core::Geometry::AbstractVolume* volume, const Core::Aabb& aabb ) {
    if ( volume != nullptr && volume->isDense() ) {

        m_mesh.loadGeometry( Core::Geometry::makeSharpBox( aabb ) );

        Core::Vector3Array tex_coords;
        tex_coords.resize( 24 );
        tex_coords.getMap() <<
            // R
            Scalar( 1 ),
            Scalar( 1 ), Scalar( 0 ), Scalar( 0 ),              // Bottom
            Scalar( 1 ), Scalar( 0 ), Scalar( 0 ), Scalar( 1 ), // Top
            Scalar( 1 ), Scalar( 1 ), Scalar( 1 ), Scalar( 1 ), // Right
            Scalar( 0 ), Scalar( 0 ), Scalar( 0 ), Scalar( 0 ), // Left
            Scalar( 1 ), Scalar( 0 ), Scalar( 0 ), Scalar( 1 ), // Floor
            Scalar( 1 ), Scalar( 1 ), Scalar( 0 ), Scalar( 0 ), // Ceil
            // G
            Scalar( 0 ), Scalar( 1 ), Scalar( 1 ), Scalar( 0 ), // Bottom
            Scalar( 0 ), Scalar( 0 ), Scalar( 1 ), Scalar( 1 ), // Top
            Scalar( 1 ), Scalar( 0 ), Scalar( 0 ), Scalar( 1 ), // Right
            Scalar( 0 ), Scalar( 1 ), Scalar( 1 ), Scalar( 0 ), // Left
            Scalar( 0 ), Scalar( 0 ), Scalar( 0 ), Scalar( 0 ), // Floor
            Scalar( 1 ), Scalar( 1 ), Scalar( 1 ), Scalar( 1 ), // Ceil
            // B
            Scalar( 0 ), Scalar( 0 ), Scalar( 0 ), Scalar( 0 ), // Bottom
            Scalar( 1 ), Scalar( 1 ), Scalar( 1 ), Scalar( 1 ), // Top
            Scalar( 0 ), Scalar( 0 ), Scalar( 1 ), Scalar( 1 ), // Right
            Scalar( 0 ), Scalar( 0 ), Scalar( 1 ), Scalar( 1 ), // Left
            Scalar( 0 ), Scalar( 0 ), Scalar( 1 ), Scalar( 1 ), // Floor
            Scalar( 0 ), Scalar( 1 ), Scalar( 1 ), Scalar( 0 ); // Ceil
        m_mesh.addAttrib( Ra::Core::Geometry::getAttribName( Ra::Core::Geometry::VERTEX_TEXCOORD ),
                          tex_coords );

        Core::Geometry::VolumeGrid* grid = static_cast<Core::Geometry::VolumeGrid*>( volume );
        m_volume = std::unique_ptr<Core::Geometry::AbstractVolume>( volume );

        auto dim = grid->size();
        TextureParameters texparam { getName(),
                                     GL_TEXTURE_3D,
                                     size_t( dim( 0 ) ),
                                     size_t( dim( 1 ) ),
                                     size_t( dim( 2 ) ),
                                     GL_RED,
                                     GL_R32F,
                                     GL_SCALAR,
                                     GL_CLAMP_TO_BORDER,
                                     GL_CLAMP_TO_BORDER,
                                     GL_CLAMP_TO_BORDER,
                                     GL_LINEAR,
                                     GL_LINEAR,
                                     grid->data().data() };
        m_tex.setParameters( texparam );

        m_isDirty = true;
    }
}

void VolumeObject::loadGeometry( Core::Geometry::AbstractVolume* volume ) {
    loadGeometry( volume, volume->computeAabb() );
}

void VolumeObject::updateGL() {
    if ( m_isDirty ) {
        m_mesh.updateGL();
        GL_CHECK_ERROR;
        m_tex.initializeGL();
        GL_CHECK_ERROR;
        m_isDirty = false;
    }
}

void VolumeObject::render( const ShaderProgram* prog ) {
    GL_CHECK_ERROR;
    // Cull faces

    GLboolean cullEnable = glIsEnabled( GL_CULL_FACE );
    int culledFaces;
    glGetIntegerv( GL_CULL_FACE_MODE, &culledFaces );
    int frontFaces;
    glGetIntegerv( GL_FRONT_FACE, &frontFaces );
    glCullFace( GL_BACK );
    glEnable( GL_CULL_FACE );

    m_mesh.render( prog );

    glCullFace( gl::GLenum( culledFaces ) );
    glFrontFace( gl::GLenum( frontFaces ) );
    if ( !cullEnable ) glDisable( GL_CULL_FACE );
    GL_CHECK_ERROR;
}

} // namespace Data
} // namespace Engine
} // namespace Ra
