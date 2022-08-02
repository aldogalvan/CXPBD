//
// Created by agalvan-admin on 7/26/22.
//

//------------------------------------------------------------------------------
#include "chai3d.h"
#include <thread>
//------------------------------------------------------------------------------
#include <GLFW/glfw3.h>
#include "world/CXPBDDeformableObject.h"
#include "collision/CXPBDAABB.h"
#include "world/CXPBDToolMesh.h"
#include "world/CXPBDTool.h"
#include "collision/CXPBDContinuousCollisionDetection.h"
#include "tetgen.h"
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

// A deformable object using the XPBD library
cXPBDDeformableMesh* box;

// a world that contains all objects of the virtual environment
cWorld* world;

// a camera to render the world in the window display
cCamera* camera;

// a light source to illuminate the objects in the world
cDirectionalLight* light;

cPositionalLight* light2;

// a colored background
cBackground* background;

cXPBDDeformableMesh* xpbd_mesh;

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

// the lines used for gripping
cShapeLine* tool1;
cShapeLine* tool2;

// the length of the tool used for visualizetion
double toolLength = 0.1;

// this is the proxy stiffness
double proxy_stiffness = 500;

// the external force applied to the object based on collision
Eigen::MatrixXd externalForce;

// force
cVector3d force(0,0,0);

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

// this function creates the tetrahedral mesh
void createTetrahedralMesh(void);

// this function progresses time
void timestep(Eigen::MatrixXd& fext, double& timestep, std::uint32_t iterations, bool gravityEnabled);

// this function computes the collision constraints
void proxyCollision(Eigen::Vector3d& goal_ , Eigen::Vector3d& proxy_, Eigen::MatrixXd& p_, Eigen::MatrixXd& plast_,
                    cXPBDDeformableMesh* a_mesh, double t_, const ColInfo& collisions);

void implicitCollision(Eigen::Vector3d& goal_ , Eigen::Vector3d& proxy_, Eigen::MatrixXd& p_, Eigen::MatrixXd& plast_,
                       cXPBDDeformableMesh* a_mesh, double t_, const ColInfo& collisions);

void testFriction(Eigen::Vector3d& goal_ , Eigen::Vector3d& proxy_, Eigen::MatrixXd& p_, Eigen::MatrixXd& plast_,
                  cXPBDDeformableMesh* a_mesh, double t_, const ColInfo& collisions);

void implicitCollision2(Eigen::Vector3d pos_ , Eigen::MatrixXd& p_, cXPBDDeformableMesh* a_mesh,  set<int> collisions);

void drawSuture(int n);


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
    camera->set(cVector3d(0.5, 0.0, 0.5),    // camera position (eye)
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

    // creates the deformable objects
    xpbd_mesh = new cXPBDDeformableMesh();
    world->addChild(xpbd_mesh);
    createTetrahedralMesh();
    xpbd_mesh->setLocalPos(Eigen::Vector3d(0,0,-0.5));
    xpbd_mesh->scaleObject(0.2);
    xpbd_mesh->connectToChai3d();
    Eigen::MatrixXd vel(xpbd_mesh->numVerts(),3);
    vel.setZero();
    xpbd_mesh->setVelocities(vel);

    // apply edge length constraint
    xpbd_mesh->constrain_edge_lengths(0.05,0.00);

    // apply tetrahedron volume constraint
    xpbd_mesh->constrain_tetrahedron_volumes(0.0,0.00);

    // wireframe vis
    xpbd_mesh->setWireMode(true,true);

    // define the radius of the tool (sphere)
    toolRadius = 0.01;

    // add the line to the world
    world->addChild(tool1);

    // add the line to the world
    world->addChild(tool2);

    // set last positions
    xpbd_mesh->setVerticesLast();

    //finds the indices at the bottm
    auto pos = xpbd_mesh->positions();

    vector<int> indices(xpbd_mesh->numVerts());

    for (int i = 0; i < xpbd_mesh->numVerts() ; i++)
    {
        double height = pos(i,2);

        if (height < -.149)
        {
            indices.emplace_back(i);
        }
    }

    // Sets indices as fixed
    xpbd_mesh->constrain_nodes_positions(indices);

    // builds a boundary box
    xpbd_mesh->buildAABBBoundaryBox();

    // computes normals
    xpbd_mesh->computeNormals();

    // resize the external force vector
    externalForce.resize(xpbd_mesh->numVerts(),3);

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

    double dt = 0.001;

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
    cVector3d proxy;
    cMatrix3d rot;
    cVector3d force;



    // get initial position
    hapticDevice->getPosition(pos);

    // get the rotation
    hapticDevice->getRotation(rot);

    // create the line representing tool
    tool = new cShapeLine(pos , pos + toolLength*cVector3d(0,0,1));
    world->addChild(tool);

    // set color at each point
    tool->m_colorPointA.setWhite();
    tool->m_colorPointB.setWhite();

    // stiffess constant
    double k = 400;
    double b = 1;

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

        // gets the current position
        hapticDevice->getPosition(pos);

        // get the rotation
        hapticDevice->getRotation(rot);

        //std::cout << "H" << std::endl;

        // change to eigen
        Eigen::Vector3d posEigen = pos.eigen(); Eigen::Vector3d proxyEigen = proxy.eigen();

        // collision info structure
        std::vector<ColInfo*> collisions;

        // computes the proxy
        if (findCollisions(posEigen, proxyEigen, toolRadius, xpbd_mesh, collisions))
        {

            std::cout << collisions.size() << std::endl;

            for (int i = 0u; i  < collisions.size() ; i++)
            {
                Eigen::Vector3i idx = collisions[i]->triangle;
                Eigen::Vector3d force_eigen = k*(proxyEigen - posEigen);
                externalForce.row(idx(0)) += -force_eigen / 3;
                externalForce.row(idx(1)) += -force_eigen / 3;
                externalForce.row(idx(2)) += -force_eigen / 3;
                //force = force_eigen;

            }

            // update the dynamics
            timestep(externalForce, dt,5,true);


        }
        else
        {
            proxy = pos;
        }

        // sets the force equal zero
        hapticDevice->setForceAndTorqueAndGripperForce(force,cVector3d(0,0,0),0);

        // draw the tool
        tool->m_pointA = cVector3d(proxyEigen);
        tool->m_pointB = cVector3d(proxyEigen) + toolLength*cMul(rot,cVector3d(0,0,1));


        // signal frequency counter
        freqCounterHaptics.signal(1);
    }

    // exit haptics thread
    simulationFinished = true;
}

// -----------------------------------------------------------------------------

void createTetrahedralMesh(void)
{

    tetgenio input;

    // TetGen switches
    char TETGEN_SWITCHES[] = "pq1.414a0.002";

    if (input.load_off(RESOURCE_PATH("../../resources/ducky/cube.off")))
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

            points.row(p) = Eigen::RowVector3d(output.pointlist[pi + 0],
                                               output.pointlist[pi + 1],
                                               output.pointlist[pi + 2]);
        }

        // sets the vertices of the mesh
        xpbd_mesh->setVertices(points);

        Eigen::MatrixXi faces(output.numberoftrifaces,3);
        auto adjtet = output.adjtetlist;

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
        xpbd_mesh->setFaces(faces);

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

        xpbd_mesh->setInsideOutside(inside,outside);

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

        xpbd_mesh->setTetrahedra(tetrahedra);

        Eigen::VectorXi facemap(faces.rows());
        for (int i = 0; i < faces.rows() ; i++)
        {
            Eigen::RowVector3i face = faces.row(i);
            std::set<int> faceset {face(0),face(1),face(2)};
            for (int j = 0 ; j < tetrahedra.rows() ; j++)
            {
                Eigen::RowVector4i tet = tetrahedra.row(j);
                std::set<int> tetset {tet(0),tet(1),tet(2),tet(3)};

                if(std::includes(tetset.begin(), tetset.end(),faceset.begin(),faceset.end()))
                {
                    facemap(i) = j;
                }
            }
        }

        xpbd_mesh->setfacemap(facemap);

        // get all the edges of our tetrahedra
        set<pair<int, int>> springs;

        for (int t = 0, ti = 0; t < output.numberoftetrahedra; ++t, ti += 4)
        {
            // store each edge of the tetrahedron as a pair of indices
            for (int i = 0; i < 4; ++i) {
                int v0 = output.tetrahedronlist[ti + i];
                for (int j = i + 1; j < 4; ++j)
                {
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
        xpbd_mesh->setEdges(edges);
    }

    Eigen::VectorXd mass;
    mass.setOnes(xpbd_mesh->numVerts());
    mass *= .0001;
    xpbd_mesh->setMass(mass);
}

void timestep(
        Eigen::MatrixXd& fext,
        double& timestep,
        std::uint32_t iterations,
        bool gravityEnabled)
{

    // start of timestep clock
    cPrecisionClock clock1;
    clock1.start(true);

    // all object constraints
    auto const& constraints   = xpbd_mesh->constraints();

    // number of constraints
    auto const J = constraints.size();

    // vector of lagrange multipliers
    std::vector<double> lagrange_multipliers(J, 0.);

    // gets the object velocity and positions
    auto& v = xpbd_mesh->velocity();
    auto& x = xpbd_mesh->positions();

    // get the mass and acceleration
    auto const& m = xpbd_mesh->mass();
    Eigen::MatrixX3d a = fext.array().colwise() / m.array();

    // set the force as zero
    fext.setZero();

    // explicit euler step
    auto vexplicit =  v + timestep * a;
    Eigen::MatrixXd p = x + timestep * vexplicit;

    // compute a new boundary box
    xpbd_mesh->buildAABBBoundaryBox(p);

    // computes new normals
    xpbd_mesh->computeNormals(p);

    // sequential gauss seidel type solve
    std::fill(lagrange_multipliers.begin(), lagrange_multipliers.end(), 0.0);
    Eigen::Vector3d F(0,0,0);

    for (auto n = 0u; n < iterations; ++n)
    {
        for (auto j = 0u; j < J; ++j)
        {
            auto const& constraint = constraints[j];
            constraint->project(p, x, m, lagrange_multipliers[j], timestep,F);
        }
    }

    // set the last positions
    xpbd_mesh->setVerticesLast();

    // update solution
    for (auto i = 0u; i < x.rows(); ++i)
    {
        v.row(i) = (p.row(i) - x.row(i)) / timestep;
        x.row(i) = p.row(i);

    }

}

void drawSuture(int n)
{

}

void drawGripper(cVector3d a_pos, double a_angle)
{

}