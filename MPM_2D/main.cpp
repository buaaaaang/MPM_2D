//
//  main.cpp
//  MPM_2D
//
//  Created by LeeSangheon on 2024/02/06.
//

#include <cassert>

#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#define MTK_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#include <Metal/Metal.hpp>
#include <AppKit/AppKit.hpp>
#include <MetalKit/MetalKit.hpp>

#include <simd/simd.h>

#include "mpm_solver.hpp"
#include "recorder.hpp"

static constexpr size_t kMaxFramesInFlight = 3;
static constexpr int screenWidth = 800;
static constexpr int screenHeight = 800;
static constexpr bool record = false;
static MTL::ClearColor backgroundColor = MTL::ClearColor::Make( 0.067, 0.176, 0.247, 1.0 );
static NS::String* title =
    NS::String::string( "simple_mpm", NS::StringEncoding::UTF8StringEncoding );

#pragma region Declarations {

class Renderer
{
public:
    Renderer( MTL::Device* pDevice );
    ~Renderer();
    void buildShaders();
    void buildBuffers();
    void draw( MTK::View* pView );

private:
    MTL::Device* _pDevice;
    MTL::CommandQueue* _pCommandQueue;
    MTL::Library* _pShaderLibrary;
    MTL::RenderPipelineState* _pPSO;
    MTL::Buffer* _pVertexDataBuffer;
    MTL::Buffer* _pInstanceDataBuffer[kMaxFramesInFlight];
    MTL::Buffer* _pIndexBuffer;
    int _frame;
    const int nParticleEdge;
    dispatch_semaphore_t _semaphore;
    static const int kMaxFramesInFlight;
    Solver mpmSolver;
    Recorder recorder;
};

class MyMTKViewDelegate : public MTK::ViewDelegate
{
    public:
        MyMTKViewDelegate( MTL::Device* pDevice );
        virtual ~MyMTKViewDelegate() override;
        virtual void drawInMTKView( MTK::View* pView ) override;

    private:
        Renderer* _pRenderer;
};

class MyAppDelegate : public NS::ApplicationDelegate
{
    public:
        ~MyAppDelegate();

        NS::Menu* createMenuBar();

        virtual void applicationWillFinishLaunching( NS::Notification* pNotification ) override;
        virtual void applicationDidFinishLaunching( NS::Notification* pNotification ) override;
        virtual bool applicationShouldTerminateAfterLastWindowClosed( NS::Application* pSender ) override;

    private:
        NS::Window* _pWindow;
        MTK::View* _pMtkView;
        MTL::Device* _pDevice;
        MyMTKViewDelegate* _pViewDelegate = nullptr;
};

#pragma endregion Declarations }


int main( int argc, char* argv[] )
{
    NS::AutoreleasePool* pAutoreleasePool = NS::AutoreleasePool::alloc()->init();

    MyAppDelegate del;

    NS::Application* pSharedApplication = NS::Application::sharedApplication();
    pSharedApplication->setDelegate( &del );
    pSharedApplication->run();

    pAutoreleasePool->release();

    return 0;
}


#pragma mark - AppDelegate
#pragma region AppDelegate {

MyAppDelegate::~MyAppDelegate()
{
    _pMtkView->release();
    _pWindow->release();
    _pDevice->release();
    delete _pViewDelegate;
}

NS::Menu* MyAppDelegate::createMenuBar()
{
    using NS::StringEncoding::UTF8StringEncoding;

    NS::Menu* pMainMenu = NS::Menu::alloc()->init();
    NS::MenuItem* pAppMenuItem = NS::MenuItem::alloc()->init();
    NS::Menu* pAppMenu = NS::Menu::alloc()->init( NS::String::string( "Appname", UTF8StringEncoding ) );

    NS::String* appName = NS::RunningApplication::currentApplication()->localizedName();
    NS::String* quitItemName = NS::String::string( "Quit ", UTF8StringEncoding )->stringByAppendingString( appName );
    SEL quitCb = NS::MenuItem::registerActionCallback( "appQuit", [](void*,SEL,const NS::Object* pSender){
        auto pApp = NS::Application::sharedApplication();
        pApp->terminate( pSender );
    } );

    NS::MenuItem* pAppQuitItem = pAppMenu->addItem( quitItemName, quitCb, NS::String::string( "q", UTF8StringEncoding ) );
    pAppQuitItem->setKeyEquivalentModifierMask( NS::EventModifierFlagCommand );
    pAppMenuItem->setSubmenu( pAppMenu );

    NS::MenuItem* pWindowMenuItem = NS::MenuItem::alloc()->init();
    NS::Menu* pWindowMenu = NS::Menu::alloc()->init( NS::String::string( "Window", UTF8StringEncoding ) );

    SEL closeWindowCb = NS::MenuItem::registerActionCallback( "windowClose", [](void*, SEL, const NS::Object*){
        auto pApp = NS::Application::sharedApplication();
            pApp->windows()->object< NS::Window >(0)->close();
    } );
    NS::MenuItem* pCloseWindowItem = pWindowMenu->addItem( NS::String::string( "Close Window", UTF8StringEncoding ), closeWindowCb, NS::String::string( "w", UTF8StringEncoding ) );
    pCloseWindowItem->setKeyEquivalentModifierMask( NS::EventModifierFlagCommand );

    pWindowMenuItem->setSubmenu( pWindowMenu );

    pMainMenu->addItem( pAppMenuItem );
    pMainMenu->addItem( pWindowMenuItem );

    pAppMenuItem->release();
    pWindowMenuItem->release();
    pAppMenu->release();
    pWindowMenu->release();

    return pMainMenu->autorelease();
}

void MyAppDelegate::applicationWillFinishLaunching( NS::Notification* pNotification )
{
    NS::Menu* pMenu = createMenuBar();
    NS::Application* pApp = reinterpret_cast< NS::Application* >( pNotification->object() );
    pApp->setMainMenu( pMenu );
    pApp->setActivationPolicy( NS::ActivationPolicy::ActivationPolicyRegular );
}

void MyAppDelegate::applicationDidFinishLaunching( NS::Notification* pNotification )
{
    CGRect frame = (CGRect){ {100.0, 100.0}, {screenWidth, screenHeight} };

    _pWindow = NS::Window::alloc()->init(
        frame,
        NS::WindowStyleMaskClosable|NS::WindowStyleMaskTitled,
        NS::BackingStoreBuffered,
        false );

    _pDevice = MTL::CreateSystemDefaultDevice();

    _pMtkView = MTK::View::alloc()->init( frame, _pDevice );
    _pMtkView->setColorPixelFormat( MTL::PixelFormat::PixelFormatBGRA8Unorm );
    _pMtkView->setClearColor( backgroundColor );

    _pViewDelegate = new MyMTKViewDelegate( _pDevice );
    _pMtkView->setDelegate( _pViewDelegate );
    
    _pMtkView->setFramebufferOnly(false); /// added to save video

    _pWindow->setContentView( _pMtkView );
    _pWindow->setTitle( title );

    _pWindow->makeKeyAndOrderFront( nullptr );

    NS::Application* pApp = reinterpret_cast< NS::Application* >( pNotification->object() );
    pApp->activateIgnoringOtherApps( true );
}

bool MyAppDelegate::applicationShouldTerminateAfterLastWindowClosed( NS::Application* pSender )
{
    return true;
}

#pragma endregion AppDelegate }


#pragma mark - ViewDelegate
#pragma region ViewDelegate {

MyMTKViewDelegate::MyMTKViewDelegate( MTL::Device* pDevice )
: MTK::ViewDelegate()
, _pRenderer( new Renderer( pDevice ) )
{
}

MyMTKViewDelegate::~MyMTKViewDelegate()
{
    delete _pRenderer;
}

void MyMTKViewDelegate::drawInMTKView( MTK::View* pView )
{
    _pRenderer->draw( pView );
}

#pragma endregion ViewDelegate }


#pragma mark - Renderer
#pragma region Renderer {

const int Renderer::kMaxFramesInFlight = 3;

Renderer::Renderer( MTL::Device* pDevice )
: _pDevice( pDevice->retain() ), _frame( 0 ), nParticleEdge(10),
    mpmSolver(), recorder(screenWidth, screenHeight)
{
    _pCommandQueue = _pDevice->newCommandQueue();
    buildShaders();
    buildBuffers();

    _semaphore = dispatch_semaphore_create( Renderer::kMaxFramesInFlight );
}

Renderer::~Renderer()
{
    _pShaderLibrary->release();
    _pVertexDataBuffer->release();
    for ( int i = 0; i < kMaxFramesInFlight; ++i )
    {
        _pInstanceDataBuffer[i]->release();
    }
    _pIndexBuffer->release();
    _pPSO->release();
    _pCommandQueue->release();
    _pDevice->release();
}

namespace shader_types
{
    struct InstanceData
    {
        simd::float4x4 instanceTransform;
        simd::float4 instanceColor;
    };
}

void Renderer::buildShaders()
{
    using NS::StringEncoding::UTF8StringEncoding;

    const char* shaderSrc = R"(
        #include <metal_stdlib>
        using namespace metal;

        struct v2f
        {
            float4 position [[position]];
            half3 color;
        };

        struct VertexData
        {
            float3 position;
        };

        struct InstanceData
        {
            float4x4 instanceTransform;
            float4 instanceColor;
        };

        v2f vertex vertexMain( device const VertexData* vertexData [[buffer(0)]],
                               device const InstanceData* instanceData [[buffer(1)]],
                               uint vertexId [[vertex_id]],
                               uint instanceId [[instance_id]] )
        {
            v2f o;
            float4 pos = float4( vertexData[ vertexId ].position, 1.0 );
            o.position = instanceData[ instanceId ].instanceTransform * pos;
            o.color = half3( instanceData[ instanceId ].instanceColor.rgb );
            return o;
        }

        half4 fragment fragmentMain( v2f in [[stage_in]] )
        {
            return half4( in.color, 1.0 );
        }
    )";

    NS::Error* pError = nullptr;
    MTL::Library* pLibrary = _pDevice->newDefaultLibrary();
    
    if ( !pLibrary )
    {
        __builtin_printf( "%s", pError->localizedDescription()->utf8String() );
        assert( false );
    }

    MTL::Function* pVertexFn = pLibrary->newFunction( NS::String::string("vertexMain", UTF8StringEncoding) );
    MTL::Function* pFragFn = pLibrary->newFunction( NS::String::string("fragmentMain", UTF8StringEncoding) );

    MTL::RenderPipelineDescriptor* pDesc = MTL::RenderPipelineDescriptor::alloc()->init();
    pDesc->setVertexFunction( pVertexFn );
    pDesc->setFragmentFunction( pFragFn );
    pDesc->colorAttachments()->object(0)->setPixelFormat( MTL::PixelFormat::PixelFormatBGRA8Unorm );

    _pPSO = _pDevice->newRenderPipelineState( pDesc, &pError );
    if ( !_pPSO )
    {
        __builtin_printf( "%s", pError->localizedDescription()->utf8String() );
        assert( false );
    }

    pVertexFn->release();
    pFragFn->release();
    pDesc->release();
    _pShaderLibrary = pLibrary;
}

void Renderer::buildBuffers()
{
    using simd::float3;

    const float s = 0.5f;

//    float3 verts[] = {
//        { -s, -s, +s },
//        { +s, -s, +s },
//        { +s, +s, +s },
//        { -s, +s, +s }
//    };
//
//    uint16_t indices[] = {
//        0, 1, 2,
//        2, 3, 0,
//    };
    
    float3 verts[nParticleEdge+1];
    for (int i=0; i<nParticleEdge; i++) {
        verts[i][0] = s*sin(i*2*M_PI/nParticleEdge);
        verts[i][1] = s*cos(i*2*M_PI/nParticleEdge);
        verts[i][2] = 0;
    }
    verts[nParticleEdge][0] = 0; verts[nParticleEdge][1] = 0; verts[nParticleEdge][2] = 0;
    
    uint16_t indices[nParticleEdge*3];
    for (int i=0; i<nParticleEdge; i++) {
        indices[i*3+0] = nParticleEdge;
        indices[i*3+1] = i;
        indices[i*3+2] = (i+1)%nParticleEdge;
    }
    
    const size_t vertexDataSize = sizeof( verts );
    const size_t indexDataSize = sizeof( indices );

    MTL::Buffer* pVertexBuffer = _pDevice->newBuffer( vertexDataSize, MTL::ResourceStorageModeManaged );
    MTL::Buffer* pIndexBuffer = _pDevice->newBuffer( indexDataSize, MTL::ResourceStorageModeManaged );

    _pVertexDataBuffer = pVertexBuffer;
    _pIndexBuffer = pIndexBuffer;

    memcpy( _pVertexDataBuffer->contents(), verts, vertexDataSize );
    memcpy( _pIndexBuffer->contents(), indices, indexDataSize );

    _pVertexDataBuffer->didModifyRange( NS::Range::Make( 0, _pVertexDataBuffer->length() ) );
    _pIndexBuffer->didModifyRange( NS::Range::Make( 0, _pIndexBuffer->length() ) );

    const size_t instanceDataSize =
        kMaxFramesInFlight * mpmSolver.numParticles * sizeof( shader_types::InstanceData );
    for ( size_t i = 0; i < kMaxFramesInFlight; ++i )
    {
        _pInstanceDataBuffer[ i ] = _pDevice->newBuffer( instanceDataSize, MTL::ResourceStorageModeManaged );
    }
}

void Renderer::draw( MTK::View* pView )
{
    using simd::float4;
    using simd::float4x4;

    NS::AutoreleasePool* pPool = NS::AutoreleasePool::alloc()->init();

    _frame = (_frame + 1) % Renderer::kMaxFramesInFlight;
    MTL::Buffer* pInstanceDataBuffer = _pInstanceDataBuffer[ _frame ];

    MTL::CommandBuffer* pCmd = _pCommandQueue->commandBuffer();
    dispatch_semaphore_wait( _semaphore, DISPATCH_TIME_FOREVER );
    Renderer* pRenderer = this;
    pCmd->addCompletedHandler( ^void( MTL::CommandBuffer* pCmd ){
        dispatch_semaphore_signal( pRenderer->_semaphore );
    });
    
    mpmSolver.nextFrame();

    const float scl = 0.01f;
    shader_types::InstanceData* pInstanceData = reinterpret_cast< shader_types::InstanceData *>( pInstanceDataBuffer->contents() );
    for ( int i = 0; i < mpmSolver.numParticles; i++ )
    {
        float x = mpmSolver.getParticleX(i);
        float y = mpmSolver.getParticleY(i);
        pInstanceData[ i ].instanceTransform = (float4x4){ (float4){ scl, 0.f, 0.f, 0.f },
                                                           (float4){ 0.f, scl, 0.f, 0.f },
                                                           (float4){ 0.f, 0.f, 0.f, 0.f },
                                                           (float4){ 2*x-1, 2*y-1, 0.f, 1.f } };

        pInstanceData[ i ].instanceColor = mpmSolver.getParticleColor(i);
    }
    pInstanceDataBuffer->didModifyRange( NS::Range::Make( 0, pInstanceDataBuffer->length() ) );


    MTL::RenderPassDescriptor* pRpd = pView->currentRenderPassDescriptor();
    MTL::RenderCommandEncoder* pEnc = pCmd->renderCommandEncoder( pRpd );

    pEnc->setRenderPipelineState( _pPSO );
    pEnc->setVertexBuffer( _pVertexDataBuffer, /* offset */ 0, /* index */ 0 );
    pEnc->setVertexBuffer( pInstanceDataBuffer, /* offset */ 0, /* index */ 1 );

    pEnc->drawIndexedPrimitives( MTL::PrimitiveType::PrimitiveTypeTriangle,
                                3*nParticleEdge, MTL::IndexType::IndexTypeUInt16,
                                _pIndexBuffer,
                                0,
                                mpmSolver.numParticles );

    pEnc->endEncoding();
    pCmd->presentDrawable( pView->currentDrawable() );
    pCmd->commit();
    pPool->release();
    
    if (record) recorder.addFrame(pView->currentDrawable()->texture(), mpmSolver.frame_time);
    
}

#pragma endregion Renderer }
