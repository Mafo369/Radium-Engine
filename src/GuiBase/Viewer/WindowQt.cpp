#include "WindowQt.hpp"

#include <QApplication>
#include <QDebug>
#include <QOpenGLContext>
#include <QResizeEvent>
#include <QSurfaceFormat>

#include <Core/Utils/Log.hpp>

using namespace Ra::Core::Utils; // log

namespace Ra {
namespace Gui {

WindowQt* WindowQt::s_getProcAddressHelper = nullptr;

QSurfaceFormat defaultFormat() {
    QSurfaceFormat format;
    format.setProfile( QSurfaceFormat::CoreProfile );
#ifndef NDEBUG
    format.setOption( QSurfaceFormat::DebugContext );
#endif
    return format;
}

WindowQt::WindowQt( QScreen* screen ) :
    QWindow( screen ),
    m_context( nullptr ),
    m_updatePending( false ),
    m_glInitialized( false ) {

    m_context.reset( new QOpenGLContext() );

    if ( !s_getProcAddressHelper )
    {
        s_getProcAddressHelper = this;
    }

    // Surface format set in BaseApplication

    setSurfaceType( OpenGLSurface );
    create();

    m_context->setFormat( QSurfaceFormat::defaultFormat() );
    if ( !m_context->create() )
    {
        LOG( logINFO ) << "Could not create OpenGL context.";
        QApplication::quit();
    }

    // cleanup connection is set in BaseApplication
}

WindowQt::~WindowQt() {
    // cannot deinitialize OpenGL here as it would require the call of a virtual member function
}

QOpenGLContext* WindowQt::context() {
    return m_context.get();
}

void WindowQt::makeCurrent() {
    m_context->makeCurrent( this );
}

void WindowQt::doneCurrent() {
    m_context->doneCurrent();
}

void WindowQt::resizeEvent( QResizeEvent* event ) {
    resize( event );
    initialize();
}

void WindowQt::exposeEvent( QExposeEvent* ) {
    initialize();
}

void WindowQt::initialize() {
    if ( !m_glInitialized.load() )
    {
        makeCurrent();

        m_glInitialized = initializeGL();

        doneCurrent();
    }
}

void WindowQt::resize( QResizeEvent* event ) {
    initialize();

    makeCurrent();

    QResizeEvent deviceSpecificResizeEvent( event->size() * devicePixelRatio(),
                                            event->oldSize() * devicePixelRatio() );

    resizeGL( &deviceSpecificResizeEvent );

    doneCurrent();
}
/// paint is done by main rendering loop, initialize instead

bool WindowQt::event( QEvent* event ) {
    switch ( event->type() )
    {
    case QEvent::UpdateRequest:
        //        paint();
        return true;

    case QEvent::Enter:
        enterEvent( event );
        return true;

    case QEvent::Leave:
        leaveEvent( event );
        return true;

    default:
        return QWindow::event( event );
    }
}

bool WindowQt::initializeGL() {
    return false;
}

void WindowQt::cleanupGL() {
    if ( m_glInitialized.load() )
    {
        makeCurrent();

        deinitializeGL();

        doneCurrent();
    }
}

void WindowQt::deinitializeGL() {}

void WindowQt::resizeGL( QResizeEvent* ) {}

/// paintgl is done by main rendering loop

void WindowQt::enterEvent( QEvent* ) {}

void WindowQt::leaveEvent( QEvent* ) {}

glbinding::ProcAddress WindowQt::getProcAddress( const char* name ) {
    if ( !s_getProcAddressHelper || name == nullptr )
    {
        return nullptr;
    }

    const auto symbol = std::string( name );

#if ( QT_VERSION >= QT_VERSION_CHECK( 5, 4, 0 ) )
    const auto qtSymbol = QByteArray::fromStdString( symbol );
#else
    const auto qtSymbol = QByteArray::fromRawData( symbol.c_str(), symbol.size() );
#endif
    return s_getProcAddressHelper->m_context->getProcAddress( qtSymbol );
}

} // namespace Gui
} // namespace Ra
