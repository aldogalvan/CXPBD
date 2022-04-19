//------------------------------------------------------------------------------
#include "chai3d.h"
//------------------------------------------------------------------------------
#include "../extras/GLFW/include/GLFW/glfw3.h"
#include "../../src/CXPBDDeformableObject.h"
#include "../../src/CXPBDSolver.h"
#include "../../src/CXPBDTool.h"
#include "../../external/libigl-main/include/igl/readOBJ.h"
#include "tetgen.h"
#include <Eigen/Dense>
#include <set>

//------------------------------------------------------------------------------
using namespace chai3d;
using namespace std;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// GENERAL SETTINGS
//------------------------------------------------------------------------------

// stereo Mode
/*
    C_STEREO_DISABLED:            Stereo is disabled
    C_STEREO_ACTIVE:              Active stereo for OpenGL NVDIA QUADRO cards
    C_STEREO_PASSIVE_LEFT_RIGHT:  Passive stereo where L/R images are rendered next to each other
    C_STEREO_PASSIVE_TOP_BOTTOM:  Passive stereo where L/R images are rendered above each other
*/
cStereoMode stereoMode = C_STEREO_DISABLED;

// fullscreen mode
bool fullscreen = false;

// mirrored display
bool mirroredDisplay = false;


//------------------------------------------------------------------------------
// DECLARED VARIABLES
//------------------------------------------------------------------------------

// XPBD

// A deformable object using the XPBD library
cXPBDDeformableMesh* xpbd_mesh;

// a virtual tool representing the haptic device in the scene
cXPBDTool* tool;

// a world that contains all objects of the virtual environment
cWorld* world;

// a camera to render the world in the window display
cCamera* camera;

// a light source to illuminate the objects in the world
cSpotLight *light;

// a colored background
cBackground* background;

// a font for rendering text
cFontPtr font;

// a label to display the rate [Hz] at which the simulation is running
cLabel* labelRates;

// a flag that indicates if the haptic simulation is currently running
bool simulationRunning = false;

// a flag that indicates if the haptic simulation has terminated
bool simulationFinished = true;

// a frequency counter to measure the simulation graphic rate
cFrequencyCounter freqCounterGraphics;

// a frequency counter to measure the simulation haptic rate
cFrequencyCounter freqCounterHaptics;

// haptic thread
cThread* hapticsThread;

// Graphics thread
cThread* graphicsThread;

// a handle to window display context
GLFWwindow* window = NULL;

// current width of window
int width = 0;

// current height of window
int height = 0;

// swap interval for the display context (vertical synchronization)
int swapInterval = 1;

// root resource path
string resourceRoot;

// a haptic device handler
cHapticDeviceHandler* handler;

// a pointer to the current haptic device
cGenericHapticDevicePtr hapticDevice;

// a generic tool
cGenericTool* haptic_tool;

//------------------------------------------------------------------------------
// DECLARED FUNCTIONS
//------------------------------------------------------------------------------

// function to create the tetrahedral mesh
void createTetrahedralMesh(cXPBDDeformableMesh* a_xpbdMesh);

// function to create the tool
void importSphere(cXPBDTool* a_tool);

// callback when the window display is resized
void windowSizeCallback(GLFWwindow* a_window, int a_width, int a_height);

// callback when an error GLFW occurs
void errorCallback(int error, const char* a_description);

// callback when a key is pressed
void keyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods);

// this function renders the scene
void updateGraphics(void);

// this function contains the main haptics simulation loop
void updateHaptics(void);

// this function closes the application
void close(void);

//---------------------------------------------------------------------------
// DECLARED MACROS
//---------------------------------------------------------------------------

// convert to resource path
#define RESOURCE_PATH(p)    (char*)((resourceRoot+string(p)).c_str())

int main(int argc, char* argv[])
{
    //--------------------------------------------------------------------------
    // INITIALIZATION
    //--------------------------------------------------------------------------

    cout << endl;
    cout << "-----------------------------------" << endl;
    cout << "CHAI3D" << endl;
    cout << "Demo: 13-primitives" << endl;
    cout << "Copyright 2003-2016" << endl;
    cout << "-----------------------------------" << endl << endl << endl;
    cout << "Keyboard Options:" << endl << endl;
    cout << "[s] - Save copy of shadowmap to file" << endl;
    cout << "[f] - Enable/Disable full screen mode" << endl;
    cout << "[m] - Enable/Disable vertical mirroring" << endl;
    cout << "[q] - Exit application" << endl;
    cout << endl << endl;


    //--------------------------------------------------------------------------
    // OPEN GL - WINDOW DISPLAY
    //--------------------------------------------------------------------------

    // initialize GLFW library
    if (!glfwInit())
    {
        cout << "failed initialization" << endl;
        cSleepMs(1000);
        return 1;
    }

    // set error callback
    glfwSetErrorCallback(errorCallback);

    // compute desired size of window
    const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    int w = 0.8 * mode->height;
    int h = 0.5 * mode->height;
    int x = 0.5 * (mode->width - w);
    int y = 0.5 * (mode->height - h);

    // set OpenGL version
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

    // set active stereo mode
    if (stereoMode == C_STEREO_ACTIVE)
    {
        glfwWindowHint(GLFW_STEREO, GL_TRUE);
    }
    else
    {
        glfwWindowHint(GLFW_STEREO, GL_FALSE);
    }

    // create display context
    window = glfwCreateWindow(w, h, "CHAI3D", NULL, NULL);
    if (!window)
    {
        cout << "failed to create window" << endl;
        cSleepMs(1000);
        glfwTerminate();
        return 1;
    }

    // get width and height of window
    glfwGetWindowSize(window, &width, &height);

    // set position of window
    glfwSetWindowPos(window, x, y);

    // set key callback
    glfwSetKeyCallback(window, keyCallback);

    // set resize callback
    glfwSetWindowSizeCallback(window, windowSizeCallback);

    // set current display context
    glfwMakeContextCurrent(window);

    // sets the swap interval for the current display context
    glfwSwapInterval(swapInterval);

#ifdef GLEW_VERSION
    // initialize GLEW library
    if (glewInit() != GLEW_OK)
    {
        cout << "failed to initialize GLEW library" << endl;
        glfwTerminate();
        return 1;
    }
#endif


    //--------------------------------------------------------------------------
    // WORLD - CAMERA - LIGHTING
    //--------------------------------------------------------------------------

    // create a new world.
    world = new cWorld();

    // set the background color of the environment
    world->m_backgroundColor.setWhite();

    // create a camera and insert it into the virtual world
    camera = new cCamera(world);
    world->addChild(camera);

    // position and orient the camera
    camera->set(cVector3d(0.9, 0.0, 0.6),    // camera position (eye)
                cVector3d(0.0, 0.0, 0.0),    // lookat position (target)
                cVector3d(0.0, 0.0, 1.0));   // direction of the (up) vector

    // set the near and far clipping planes of the camera
    // anything in front or behind these clipping planes will not be rendered
    camera->setClippingPlanes(0.01, 10.0);

    // set stereo mode
    camera->setStereoMode(stereoMode);

    // set stereo eye separation and focal length (applies only if stereo is enabled)
    camera->setStereoEyeSeparation(0.03);
    camera->setStereoFocalLength(1.8);

    // set vertical mirrored display mode
    camera->setMirrorVertical(mirroredDisplay);

    // create a light source
    light = new cSpotLight(world);

    // attach light to camera
    world->addChild(light);

    // enable light source
    light->setEnabled(true);

    // position the light source
    light->setLocalPos(0.6, 0.6, 0.5);

    // define the direction of the light beam
    light->setDir(-0.5,-0.5,-0.5);

    // enable this light source to generate shadows
    light->setShadowMapEnabled(true);

    // set the resolution of the shadow map
    //light->m_shadowMap->setQualityLow();
    light->m_shadowMap->setQualityMedium();

    // set light cone half angle
    light->setCutOffAngleDeg(30);

    //--------------------------------------------------------------------------
    // HAPTIC DEVICES / TOOLS
    //--------------------------------------------------------------------------

    // create a haptic device handler
    handler = new cHapticDeviceHandler();

    // get access to the first available haptic device found
    handler->getDevice(hapticDevice, 0);

    // retrieve information about the current haptic device
    cHapticDeviceInfo hapticDeviceInfo = hapticDevice->getSpecifications();



    //--------------------------------------------------------------------------
    // CREATE XPBD OBJECT
    //--------------------------------------------------------------------------

    // creates the deformable objects
    xpbd_mesh = new cXPBDDeformableMesh();
    createTetrahedralMesh(xpbd_mesh);
    xpbd_mesh->setLocalPos(Eigen::Vector3d(0,0,5));
    xpbd_mesh->scaleObject(0.1);
    xpbd_mesh->connectToChai3d();

    Eigen::MatrixXd vel(xpbd_mesh->numVerts(),3);
    vel.setZero();
    xpbd_mesh->setVelocities(vel);
    //chai3d_mesh->computeAllEdges();
    //chai3d_mesh->setShowEdges(true);
    //chai3d_mesh->m_edgeLineColor = cColorf(1,0,0,0);
    //chai3d_mesh->setShowTriangles(true);

    xpbd_mesh->constrain_edge_lengths(0.01);
    xpbd_mesh->constrain_tetrahedron_volumes(0.01);
    //xpbd_mesh->constrain_neohookean_elasticity_potential(100000,0.5);

    // creates the tool
    tool = new cXPBDTool(world);
    // connect the haptic device to the virtual tool
    tool->m_tool->setHapticDevice(hapticDevice);

    importSphere(tool);
    world->addChild(tool);
    world->addChild(xpbd_mesh);


    // set last positions
    xpbd_mesh->setVerticesLast();
    tool->setVerticesLast();

    //--------------------------------------------------------------------------
    // WIDGETS
    //--------------------------------------------------------------------------

    // create a font
    font = NEW_CFONTCALIBRI20();

    // create a label to display the haptic and graphic rate of the simulation
    labelRates = new cLabel(font);
    labelRates->m_fontColor.setBlack();
    camera->m_frontLayer->addChild(labelRates);

    // create a background
    background = new cBackground();
    camera->m_backLayer->addChild(background);

    // set background properties
    background->setCornerColors(cColorf(1.0f, 1.0f, 1.0f),
                                cColorf(1.0f, 1.0f, 1.0f),
                                cColorf(0.8f, 0.8f, 0.8f),
                                cColorf(0.8f, 0.8f, 0.8f));


    //--------------------------------------------------------------------------
    // START SIMULATION
    //--------------------------------------------------------------------------

    // create a thread which starts the main haptics rendering loop
    hapticsThread = new cThread();
    hapticsThread->start(updateHaptics, CTHREAD_PRIORITY_HAPTICS);

    // creates a thread which starts the main graphics rendering loop
    //graphicsThread = new cThread();
    //graphicsThread->start(updateGraphics, CTHREAD_PRIORITY_GRAPHICS);

    // setup callback when application exits
    atexit(close);

    //--------------------------------------------------------------------------
    // MAIN GRAPHIC LOOP
    //--------------------------------------------------------------------------

    // call window size callback at initialization
    windowSizeCallback(window, width, height);

    // main graphic loop
    while (!glfwWindowShouldClose(window))
    {
        // get width and height of window
        glfwGetWindowSize(window, &width, &height);

        // Solver
        solve(xpbd_mesh, tool, Eigen::MatrixXd::Zero(xpbd_mesh->numVerts(),
                                                     3),0.01,10,5,true);

        updateGraphics();

        // swap buffers
        glfwSwapBuffers(window);

        // process events
        glfwPollEvents();

        // signal frequency counter
        freqCounterGraphics.signal(1);


    }

    // close window
    glfwDestroyWindow(window);

    // terminate GLFW library
    glfwTerminate();

    // exit
    return (0);
}

//------------------------------------------------------------------------------

void windowSizeCallback(GLFWwindow* a_window, int a_width, int a_height)
{
    // update window size
    width = a_width;
    height = a_height;
}

//------------------------------------------------------------------------------

void errorCallback(int a_error, const char* a_description)
{
    cout << "Error: " << a_description << endl;
}

//------------------------------------------------------------------------------

void keyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods)
{
    // filter calls that only include a key press
    if ((a_action != GLFW_PRESS) && (a_action != GLFW_REPEAT))
    {
        return;
    }

        // option - exit
    else if ((a_key == GLFW_KEY_ESCAPE) || (a_key == GLFW_KEY_Q))
    {
        glfwSetWindowShouldClose(a_window, GLFW_TRUE);
    }

        // option - save shadow map to file
    else if (a_key == GLFW_KEY_S)
    {
        cImagePtr image = cImage::create();
        light->m_shadowMap->copyDepthBuffer(image);
        image->saveToFile("shadowmapshot.png");
        cout << "> Saved screenshot of shadowmap to file.       \r";
    }

        // option - toggle fullscreen
    else if (a_key == GLFW_KEY_F)
    {
        // toggle state variable
        fullscreen = !fullscreen;

        // get handle to monitor
        GLFWmonitor* monitor = glfwGetPrimaryMonitor();

        // get information about monitor
        const GLFWvidmode* mode = glfwGetVideoMode(monitor);

        // set fullscreen or window mode
        if (fullscreen)
        {
            glfwSetWindowMonitor(window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
            glfwSwapInterval(swapInterval);
        }
        else
        {
            int w = 0.8 * mode->height;
            int h = 0.5 * mode->height;
            int x = 0.5 * (mode->width - w);
            int y = 0.5 * (mode->height - h);
            glfwSetWindowMonitor(window, NULL, x, y, w, h, mode->refreshRate);
            glfwSwapInterval(swapInterval);
        }
    }

        // option - toggle vertical mirroring
    else if (a_key == GLFW_KEY_M)
    {
        mirroredDisplay = !mirroredDisplay;
        camera->setMirrorVertical(mirroredDisplay);
    }
}

//------------------------------------------------------------------------------

void close(void)
{
    // stop the simulation
    simulationRunning = false;

    // wait for graphics and haptics loops to terminate
    while (!simulationFinished) { cSleepMs(100); }

    // delete resources
    delete hapticsThread;
    delete world;

}

//------------------------------------------------------------------------------

void updateGraphics(void)
{
    /////////////////////////////////////////////////////////////////////
    // UPDATE WIDGETS
    /////////////////////////////////////////////////////////////////////

    // update haptic and graphic rate data
    labelRates->setText(cStr(freqCounterGraphics.getFrequency(), 0) + " Hz / " +
                        cStr(freqCounterHaptics.getFrequency(), 0) + " Hz");

    // update position of label
    labelRates->setLocalPos((int)(0.5 * (width - labelRates->getWidth())), 15);

    /////////////////////////////////////////////////////////////////////
    // RENDER SCENE
    /////////////////////////////////////////////////////////////////////

    // update shadow maps (if any)
    world->updateShadowMaps(false, mirroredDisplay);

    // render world
    camera->renderView(width, height);

    // wait until all GL commands are completed
    glFinish();

    // check for any OpenGL errors
    GLenum err = glGetError();
    if (err != GL_NO_ERROR) cout << "Error: " << gluErrorString(err) << endl;
}

//------------------------------------------------------------------------------

enum cMode
{
    IDLE,
    SELECTION
};

void updateHaptics(void)
{
    cMode state = IDLE;
    cGenericObject* object = NULL;
    cTransform tool_T_object;

    // simulation in now running
    simulationRunning  = true;
    simulationFinished = false;

    // main haptic simulation loop
    while(simulationRunning)
    {
        /////////////////////////////////////////////////////////////////////////
        // HAPTIC RENDERING
        /////////////////////////////////////////////////////////////////////////

        // signal frequency counter
        freqCounterHaptics.signal(1);

        // compute global reference frames for each object
        world->computeGlobalPositions(true);

        // update position and orientation of tool
        tool->m_tool->updateFromDevice();

        // update the visualizer


    }

    // exit haptics thread
    simulationFinished = true;
}

void createTetrahedralMesh(cXPBDDeformableMesh* a_xpbdMesh)
{

    tetgenio input;

    // TetGen switches
    char TETGEN_SWITCHES[] = "pq1.414a0.002";

    if (input.load_off(RESOURCE_PATH("../../resources/ducky/duck-200.off")))
    {
        // use TetGen to tetrahedralize our mesh
        tetgenio output;
        tetrahedralize(TETGEN_SWITCHES, &input, &output);

        Eigen::MatrixXd points(output.numberofpoints,3);

        // create a vertex in the object for each point of the result
        for (int p = 0, pi = 0; p < output.numberofpoints; ++p, pi += 3)
        {
            cVector3d point;
            point.x(output.pointlist[pi + 0]);
            point.y(output.pointlist[pi + 1]);
            point.z(output.pointlist[pi + 2]);
            //a_chai3dMesh->newVertex(point);

            points.row(p) = Eigen::RowVector3d(output.pointlist[pi + 0],
                                                       output.pointlist[pi + 1],
                                                       output.pointlist[pi + 2]);

        }

        // sets the vertices of the mesh
        a_xpbdMesh->setVertices(points);

        Eigen::MatrixXi faces(output.numberoftrifaces,3);

        // create a triangle for each face on the surface
        for (int t = 0, ti = 0; t < output.numberoftrifaces; ++t, ti += 3)
        {
            cVector3d p[3];
            unsigned int vi[3];
            for (int i = 0; i < 3; ++i)
            {
                int tc = output.trifacelist[ti + i];
                vi[i] = tc;
                int pi = tc * 3;
                p[i].x(output.pointlist[pi + 0]);
                p[i].y(output.pointlist[pi + 1]);
                p[i].z(output.pointlist[pi + 2]);
            }
            //unsigned int index = a_object->newTriangle(p[1], p[0], p[2]);
            //a_chai3dMesh->newTriangle(vi[0], vi[1], vi[2]);
            faces.row(t) = Eigen::RowVector3i(vi[0],vi[1],vi[2]);
        }

        // sets the faces of the mesh
        a_xpbdMesh->setFaces(faces);

        // find out exactly which vertices are on the inside and outside
        set<int> inside, outside;
        for (int t = 0; t < output.numberoftrifaces * 3; ++t)
        {
            outside.insert(output.trifacelist[t]);
        }
        for (int p = 0; p < output.numberofpoints; ++p)
        {
            if (outside.find(p) == outside.end())
                inside.insert(p);
        }

        a_xpbdMesh->setInsideOutside(inside,outside);

        Eigen::MatrixXi tetrahedra(output.numberoftetrahedra,4);

        for (int t = 0, ti = 0; t < output.numberoftetrahedra; ++t, ti += 4)
        {
            int v0 = output.tetrahedronlist[ti + 0];
            int v1 = output.tetrahedronlist[ti + 1];
            int v2 = output.tetrahedronlist[ti + 2];
            int v3 = output.tetrahedronlist[ti + 3];

            Eigen::RowVector4i tetrahedron;
            tetrahedron[0] = v0;
            tetrahedron[1] = v1;
            tetrahedron[2] = v2;
            tetrahedron[3] = v3;

            tetrahedra.row(t) = (tetrahedron);
        }

        a_xpbdMesh->setTetrahedra(tetrahedra);

        // get all the edges of our tetrahedra
        set<pair<int, int>> springs;

        for (int t = 0, ti = 0; t < output.numberoftetrahedra; ++t, ti += 4)
        {
            // store each edge of the tetrahedron as a pair of indices
            for (int i = 0; i < 4; ++i) {
                int v0 = output.tetrahedronlist[ti + i];
                for (int j = i + 1; j < 4; ++j) {
                    int v1 = output.tetrahedronlist[ti + j];
                    springs.insert(pair<int, int>(min(v0, v1), max(v0, v1)));
                }
            }
        }

        Eigen::MatrixXi edges(springs.size(),2);

        int i = 0;

        for (auto sp : springs)
        {
            edges.row(i) = Eigen::RowVector2i(sp.first, sp.second);
            i++;
        }

        a_xpbdMesh->setEdges(edges);
    }

    Eigen::VectorXd mass;
    mass.setOnes(a_xpbdMesh->numVerts());
    a_xpbdMesh->setMass(mass*0.0001);
}

void importSphere(cXPBDTool* a_tool)
{
    Eigen::MatrixX3d V;
    Eigen::MatrixX3i F;
    igl::readOBJ(RESOURCE_PATH("../../resources/sphere/sphere.obj"), V, F);
    a_tool->setVertices(V);
    a_tool->setFaces(F);
    a_tool->scaleObject(0.05);
    a_tool->connectToChai3d();
}
//------------------------------------------------------------------------------
