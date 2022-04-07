#include <Headless/OpenGLContext/OpenGLContext.hpp>

#include <GLFW/glfw3.h>

#include <glbinding/AbstractFunction.h>
#include <glbinding/Binding.h>
#include <glbinding/CallbackMask.h>
#include <glbinding/FunctionCall.h>
#include <glbinding/Version.h>
#include <glbinding/glbinding.h>

#include <glbinding/gl/gl.h>

#include <glbinding-aux/ContextInfo.h>
#include <glbinding-aux/Meta.h>
#include <glbinding-aux/ValidVersions.h>
#include <glbinding-aux/types_to_string.h>

#include <globjects/globjects.h>

#include <iostream>

namespace Ra {
namespace Headless {
using namespace gl;
using namespace glbinding;

static void error( int errnum, const char* errmsg ) {
    globjects::critical() << errnum << ": " << errmsg << std::endl;
}

OpenGLContext::OpenGLContext( const std::array<int, 2>& size ) {
    // initialize openGL
    if ( glfwInit() ) {
        glfwSetErrorCallback( error );
        glfwDefaultWindowHints();
        glfwWindowHint( GLFW_VISIBLE, false );
        glfwWindowHint( GLFW_CONTEXT_VERSION_MAJOR, 4 );
        glfwWindowHint( GLFW_CONTEXT_VERSION_MINOR, 1 );
        glfwWindowHint( GLFW_OPENGL_FORWARD_COMPAT, true );
        glfwWindowHint( GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE );
        m_glfwContext =
            glfwCreateWindow( size[0], size[1], "Radium CommandLine Context", nullptr, nullptr );
    }
    if ( m_glfwContext == nullptr ) {
        std::cerr << "OpenGL context creation failed. Terminate execution." << std::endl;
        const char* description;
        int code = glfwGetError( &description );
        std::cerr << "\tError code :" << code << std::endl
                  << "\t error string : " << description << std::endl;
        glfwTerminate();
        std::exit( -1 );
    }
    else {
        // Initialize globjects (internally initializes glbinding, and registers the current
        // context)
        glfwMakeContextCurrent( m_glfwContext );
        globjects::init( []( const char* name ) { return glfwGetProcAddress( name ); } );
        glfwSetWindowUserPointer( m_glfwContext, this );
        auto resizeCB = []( GLFWwindow* window, int width, int height ) {
            auto context = static_cast<OpenGLContext*>( glfwGetWindowUserPointer( window ) );
            context->resizeFrameBuffer( width, height );
        };
        glfwSetFramebufferSizeCallback( m_glfwContext, resizeCB );
        auto keyCB = []( GLFWwindow* window, int key, int scancode, int action, int mods ) {
            auto context = static_cast<OpenGLContext*>( glfwGetWindowUserPointer( window ) );
            context->keyboardEventCalback( key, scancode, action, mods );
        };
        glfwSetKeyCallback( m_glfwContext, keyCB );

        auto mouseCB = []( GLFWwindow* window, int button, int action, int mods ) {
            auto context = static_cast<OpenGLContext*>( glfwGetWindowUserPointer( window ) );
            // see https://www.glfw.org/docs/latest/window_guide.html#window_scale
            float xscale, yscale;
            glfwGetWindowContentScale( window, &xscale, &yscale );
            double xpos, ypos;
            glfwGetCursorPos( window, &xpos, &ypos );
            context->mouseEventCalback(
                button, action, mods, int( xpos / xscale ), int( ypos / yscale ) );
        };
        glfwSetMouseButtonCallback( m_glfwContext, mouseCB );

        auto scrollCB = []( GLFWwindow* window, double xoffset, double yoffset ) {
            auto context = static_cast<OpenGLContext*>( glfwGetWindowUserPointer( window ) );
            float xscale, yscale;
            glfwGetWindowContentScale( window, &xscale, &yscale );
            context->scrollEventCalback( int( xoffset / xscale ), int( yoffset / yscale ) );
        };
        glfwSetScrollCallback( m_glfwContext, scrollCB );
    }
}
OpenGLContext::~OpenGLContext() {
    glfwTerminate();
}
void OpenGLContext::makeCurrent() const {
    if ( m_glfwContext ) { glfwMakeContextCurrent( m_glfwContext ); }
}

void OpenGLContext::doneCurrent() const {
    if ( m_glfwContext ) { glfwMakeContextCurrent( nullptr ); }
}

bool OpenGLContext::isValid() const {
    return m_glfwContext != nullptr;
}

std::string OpenGLContext::getInfo() const {
    std::stringstream infoText;
    using ContextInfo = glbinding::aux::ContextInfo;
    makeCurrent();
    infoText << "*** OffScreen OpenGL context ***" << std::endl;
    infoText << "Renderer (glbinding) : " << ContextInfo::renderer() << "\n";
    infoText << "Vendor   (glbinding) : " << ContextInfo::vendor() << "\n";
    infoText << "OpenGL   (glbinding) : " << ContextInfo::version().toString() << "\n";
    infoText << "GLSL                 : " << gl::glGetString( GL_SHADING_LANGUAGE_VERSION ) << "\n";
    doneCurrent();

    return infoText.str();
}

void OpenGLContext::show( EventMode mode, float delay ) {
    m_mode  = mode;
    m_delay = delay;
    glfwShowWindow( m_glfwContext );
}

void OpenGLContext::hide() {
    glfwHideWindow( m_glfwContext );
}

void OpenGLContext::resize( const std::array<int, 2>& size ) {
    glfwSetWindowSize( m_glfwContext, size[0], size[1] );
}

bool OpenGLContext::processEvents() {
    switch ( m_mode ) {
    case EventMode::POLL:
        glfwPollEvents();
        break;
    case EventMode::WAIT:
        glfwWaitEvents();
        break;
    case EventMode::TIMEOUT:
        glfwWaitEventsTimeout( m_delay );
        break;
    default:
        glfwPollEvents();
        break;
    }
    return true;
}

void OpenGLContext::resizeFrameBuffer( int width, int height ) {
    gl::glViewport( 0, 0, width, height );
    m_resizers.notify( width, height );
}

void OpenGLContext::keyboardEventCalback( int key, int scancode, int action, int mods ) {
    m_keyboardObservers.notify( key, scancode, action, mods );
}

void OpenGLContext::mouseEventCalback( int button, int action, int mods, int x, int y ) {
    m_mouseObservers.notify( button, action, mods, x, y );
}

void OpenGLContext::scrollEventCalback( int xoffset, int yoffset ) {
    m_scrollObservers.notify( xoffset, yoffset );
}

void OpenGLContext::renderLoop( std::function<void( float )> render ) {
    double prevFrameDate = glfwGetTime();
    double curFrameDate;

    int width, height;
    glfwGetFramebufferSize( m_glfwContext, &width, &height );
    glViewport( 0, 0, width, height );

    while ( !glfwWindowShouldClose( m_glfwContext ) ) {
        curFrameDate = glfwGetTime();
        render( curFrameDate - prevFrameDate );
        prevFrameDate = curFrameDate;

        glfwSwapBuffers( m_glfwContext );
        processEvents();
    }
}

} // namespace Headless
} // namespace Ra
