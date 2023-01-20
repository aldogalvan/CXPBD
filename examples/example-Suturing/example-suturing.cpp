#include "chai3d.h"
#include <thread>
#include "simulator.h"
#include <GLFW/glfw3.h>
#include <Eigen/Core>
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
// STATES
//------------------------------------------------------------------------------
enum MouseStates
{
    MOUSE_IDLE,
    MOUSE_MOVE_CAMERA
};

enum HapticStates
{
    HAPTIC_IDLE,
    HAPTIC_SELECTION
};



//------------------------------------------------------------------------------
// DECLARED VARIABLES
//------------------------------------------------------------------------------

// simulator
simulator* sim;

// mesh representing the object
cMesh* simvis;

// a world that contains all objects of the virtual environment
cWorld* world;

// a camera to render the world in the window display
cCamera* camera;

// a light source to illuminate the objects in the world
cDirectionalLight* light;

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

// mouse state
MouseStates mouseState = MOUSE_IDLE;

// last mouse position
double mouseX, mouseY;

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

// the radius of the tool
double toolRadius;

// the line used for visualization
cShapeLine* tool;

// the length of the tool used for visualizetion
double toolLength = 0.1;

// this is the proxy stiffness
double proxy_stiffness = 500;

// the external force applied to the object based on collision
Eigen::MatrixXd externalForce;

// force
cVector3d force(0,0,0);

// scalong
float wss = 5;

//------------------------------------------------------------------------------
// DECLARED FUNCTIONS
//------------------------------------------------------------------------------

// callback when the window display is resized
void windowSizeCallback(GLFWwindow* a_window, int a_width, int a_height);

// callback when an error GLFW occurs
void errorCallback(int error, const char* a_description);

// callback when a key is pressed
void keyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods);

// callback to handle mouse click
void mouseButtonCallback(GLFWwindow* a_window, int a_button, int a_action, int a_mods);

// callback to handle mouse motion
void mouseMotionCallback(GLFWwindow* a_window, double a_posX, double a_posY);

// callback to handle mouse scroll
void mouseScrollCallback(GLFWwindow* a_window, double a_offsetX, double a_offsetY);

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

    // set mouse position callback
    glfwSetCursorPosCallback(window, mouseMotionCallback);

    // set mouse button callback
    glfwSetMouseButtonCallback(window, mouseButtonCallback);

    // set mouse scroll callback
    glfwSetScrollCallback(window, mouseScrollCallback);

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
    world->m_backgroundColor.setBlack();

    // create a camera and insert it into the virtual world
    camera = new cCamera(world);
    world->addChild(camera);

    // position and orient the camera
    camera->set(cVector3d(1.5, 0.0, 1.5),    // camera position (eye)
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
    light = new cDirectionalLight(world);

    // attach light to camera
    world->addChild(light);

    // enable light source
    light->setEnabled(true);

    // position the light source
    light->setLocalPos(0, 0, 1);

    // define the direction of the light beam
    light->setDir(0.0,0.0,0.0);

    //--------------------------------------------------------------------------
    // HAPTIC DEVICES / TOOLS
    //--------------------------------------------------------------------------

    // create a haptic device handler
    handler = new cHapticDeviceHandler();

    // get access to the first available haptic device found
    handler->getDevice(hapticDevice, 0);

    // open a connection to haptic device
    hapticDevice->open();

    // calibrate device if necessary
    hapticDevice->calibrate();

    // retrieve information about the current haptic device
    cHapticDeviceInfo info = hapticDevice->getSpecifications();

    //--------------------------------------------------------------------------
    // CREATE XPBD OBJECT
    //--------------------------------------------------------------------------

    // creates the simulator
    sim = new simulator;

    // creates the mesh
    simvis = new cMesh();
    world->addChild(simvis);

    // imports the mesh to chai3d
    auto p = sim->positions();
    auto f = sim->triangles();
    for (int pi = 0 ; pi < p.rows();pi++)
    {
        simvis->newVertex((float)p(pi,0),(float)p(pi,1),(float)p(pi,2));
    }
    for (int fi = 0 ; fi < f.rows();fi++)
    {
        simvis->newTriangle(f(fi,0),f(fi,1),f(fi,2));
    }
    simvis->setWireMode(true,true);

    // define the radius of the tool (sphere)
    toolRadius = 0.01;

    // add the line to the world
    world->addChild(tool);

    // set the external force vector equal to zero
    externalForce.setZero();

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
    //background->setColor()
    //--------------------------------------------------------------------------
    // START SIMULATION
    //--------------------------------------------------------------------------

    // creates a thread which starts the main haptic rendering loop
    hapticsThread = new cThread();
    hapticsThread->start(updateHaptics, CTHREAD_PRIORITY_HAPTICS);

    // creates a thread which starts the main graphics rendering loop
    // TODO: FIGURE THIS OUT!
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

        // process events
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
        //light->m_shadowMap->copyDepthBuffer(image);
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

void mouseButtonCallback(GLFWwindow* a_window, int a_button, int a_action, int a_mods)
{
    if (a_button == GLFW_MOUSE_BUTTON_RIGHT && a_action == GLFW_PRESS)
    {
        // store mouse position
        glfwGetCursorPos(window, &mouseX, &mouseY);

        // update mouse state
        mouseState = MOUSE_MOVE_CAMERA;
    }

    else
    {
        // update mouse state
        mouseState = MOUSE_IDLE;
    }
}

//------------------------------------------------------------------------------

void mouseMotionCallback(GLFWwindow* a_window, double a_posX, double a_posY)
{
    if (mouseState == MOUSE_MOVE_CAMERA)
    {
        // compute mouse motion
        int dx = a_posX - mouseX;
        int dy = a_posY - mouseY;
        mouseX = a_posX;
        mouseY = a_posY;

        // compute new camera angles
        double azimuthDeg = camera->getSphericalAzimuthDeg() - 0.5 * dx;
        double polarDeg = camera->getSphericalPolarDeg() - 0.5 * dy;

        // assign new angles
        camera->setSphericalAzimuthDeg(azimuthDeg);
        camera->setSphericalPolarDeg(polarDeg);
    }
}

//------------------------------------------------------------------------------

void mouseScrollCallback(GLFWwindow* a_window, double a_offsetX, double a_offsetY)
{
    double r = camera->getSphericalRadius();
    r = cClamp(r + 0.1 * a_offsetY, 0.5, 3.0);
    camera->setSphericalRadius(r);
}

//------------------------------------------------------------------------------

void close(void)
{
    // stop the simulation
    simulationRunning = false;

    // wait for graphics and haptics loops to terminate
    while (!simulationFinished) { cSleepMs(100); }

    hapticDevice->close();

    // delete resources
    delete hapticsThread;
    delete world;
    delete handler;

}

//------------------------------------------------------------------------------

void updateGraphics(void)
{

    /////////////////////////////////////////////////////////////////////
    // UPDATE WIDGETS
    /////////////////////////////////////////////////////////////////////

    //std::cout << "G" << std::endl;

    // update haptic and graphic rate data
    labelRates->setText(cStr(freqCounterGraphics.getFrequency(), 0) + " Hz / " +
                        cStr(freqCounterHaptics.getFrequency(), 0) + " Hz");

    // update position of label
    labelRates->setLocalPos((int) (0.5 * (width - labelRates->getWidth())), 15);

    /////////////////////////////////////////////////////////////////////
    // RENDER SCENE
    /////////////////////////////////////////////////////////////////////

    auto p = sim->positions();

    //std::cout << p.row(0) << std::endl;

    // TODO: Parallelize
    for (int i = 0; i < p.rows(); i++)
    {
        simvis->m_vertices->setLocalPos(i,p(i,0),p(i,1),p(i,2));
    }

    // update shadow maps (if any)
    world->updateShadowMaps(false, mirroredDisplay);

    // render world
    camera->renderView(width, height);

    // wait until all GL commands are completed
    glFinish();

    // check for any OpenGL errors
    GLenum err;
    err = glGetError();
    if (err != GL_NO_ERROR) cout << "Error:  %s\n" << gluErrorString(err);

}

//------------------------------------------------------------------------------

enum cMode
{
    IDLE,
    SELECTION
};

void updateHaptics(void)
{

    // simulation in now running
    simulationRunning  = true;
    simulationFinished = false;

    // declare some variables
    cVector3d pos;
    cVector3d vel;
    cVector3d proxy;
    cMatrix3d theta;
    cMatrix3d omega;
    cVector3d force;
    cVector3d torque;

    // get initial position
    hapticDevice->getPosition(pos);

    // get the rotation
    hapticDevice->getRotation(theta);

    // create the line representing tool
    tool = new cShapeLine(pos , pos + toolLength*cVector3d(0,0,1));
    world->addChild(tool);

    // set color at each point
    tool->m_colorPointA.setWhite();
    tool->m_colorPointB.setWhite();

    // stiffess constant
    double k = 750;
    double b = 0;

    // friction coefficient
    double us = 0.1;
    double uk = 0.1;

    // initial step
    double dt = 0.001;

    // main haptic simulation loop
    while(simulationRunning) {


        /////////////////////////////////////////////////////////////////////
        // HAPTIC FORCE COMPUTATION
        /////////////////////////////////////////////////////////////////////

        // sets the force equal zero
        force = cVector3d(0,0,0);
        torque = cVector3d(0,0,0);

        // gets the position
        hapticDevice->getPosition(pos);
        pos = wss*pos;

        // gets the velocity
        hapticDevice->getLinearVelocity(vel);

        // get the rotation
        hapticDevice->getRotation(theta);

        // change to eigen
        float ppos[3]; float pproxy[3];
        ppos[0] = 1000*((float)pos(0) + .01); ppos[1] = 1000*((float)pos(1)+.01); ppos[2] = 1000*((float)pos(2)+0.01);
        pproxy[0] = 1000*(float)pos(0); pproxy[1] = 1000*(float)pos(1); pproxy[2] = 1000*(float)pos(2);

        // command a force and update the dynamics
        float force[3] = {0,0,0};
        sim->updateDynamics(ppos,pproxy,force, 0,0.01);

        // sets the force equal zero
        hapticDevice->setForceAndTorqueAndGripperForce(cVector3d(0,0,0),cVector3d(0,0,0),0);

        // draw the tool
        tool->m_pointA = cVector3d(pproxy[0] / 1000,pproxy[1] / 1000,pproxy[2] / 1000);
        tool->m_pointB = cVector3d(pproxy[0] / 1000,pproxy[1] / 1000,pproxy[2] / 1000)
                         + toolLength*cMul(theta,cVector3d(0,0,1));

        // signal frequency counter
        freqCounterHaptics.signal(1);
    }

    // exit haptics thread
    simulationFinished = true;
}

// -----------------------------------------------------------------------------

