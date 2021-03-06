//
// Created by aldo on 5/30/22.
//

#include "chai3d.h"
#include <Eigen/Dense>
//------------------------------------------------------------------------------
#include "../extras/GLFW/include/GLFW/glfw3.h"
#include "world/CXPBDDeformableObject.h"
#include "world/CXPBDSolver.h"
#include "collision/CXPBDAABB.h"
#include "world/CXPBDToolMesh.h"
#include "world/CXPBDTool.h"
#include "collision/CXPBDContinuousCollisionDetection.h"
#include "CXPBDVolumeConstraint.cuh"
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
// DECLARED VARIABLES
//------------------------------------------------------------------------------

// XPBD

// A deformable object using the XPBD library
cXPBDDeformableMesh* xpbd_mesh;

// a virtual tool representing the haptic device in the scene
cShapeSphere* tool;

// a new tool representing the haptic device
cShapeSphere* proxy;

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

// this is the proxy stiffness
double proxy_stiffness;


//------------------------------------------------------------------------------
// DECLARED FUNCTIONS
//------------------------------------------------------------------------------

// function to create the tetrahedral mesh
void createTetrahedralMesh(cXPBDDeformableMesh* a_xpbdMesh);

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

// this function progresses time
void timestep(cXPBDDeformableMesh* model, Eigen::MatrixX3d const& fext, double timestep = 0.01,
              std::uint32_t iterations = 10, std::uint32_t substeps   = 10, bool gravityEnabled = true);

// this function computes the collision constraints
void computeCollisionConstraintsPassive(Eigen::Vector3d goal_ , Eigen::Vector3d& proxy_, Eigen::MatrixXd& p_, Eigen::MatrixXd& plast_,
                                        cXPBDDeformableMesh* a_mesh, double t_, const ColInfo& collisions);

void computeCollisionConstraintsActive(Eigen::Vector3d goal_ , Eigen::Vector3d& proxy_, Eigen::MatrixXd& p_, Eigen::MatrixXd& plast_,
                                       cXPBDDeformableMesh* a_mesh, double t_, const ColInfo& collisions);

void testFriction(Eigen::Vector3d goal_ , Eigen::Vector3d& proxy_, Eigen::MatrixXd& p_, Eigen::MatrixXd& plast_,
                  cXPBDDeformableMesh* a_mesh, double t_, const ColInfo& collisions);

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
    camera->set(cVector3d(0.9, 0.0, -0.5),    // camera position (eye)
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

    // open a connection to haptic device
    hapticDevice->open();


    // retrieve information about the current haptic device
    cHapticDeviceInfo info = hapticDevice->getSpecifications();

    //--------------------------------------------------------------------------
    // CREATE XPBD OBJECT
    //--------------------------------------------------------------------------

    // creates the deformable objects
    xpbd_mesh = new cXPBDDeformableMesh();
    world->addChild(xpbd_mesh);
    createTetrahedralMesh(xpbd_mesh);
    xpbd_mesh->setLocalPos(Eigen::Vector3d(0,0,1));
    xpbd_mesh->scaleObject(0.2);
    xpbd_mesh->connectToChai3d();
    Eigen::MatrixXd vel(xpbd_mesh->numVerts(),3);
    vel.setZero();
    xpbd_mesh->setVelocities(vel);
    xpbd_mesh->constrain_edge_lengths(0.25);
    xpbd_mesh->constrain_tetrahedron_volumes(0.25);
    //xpbd_mesh->constrain_neohookean_elasticity_potential(10000000,0.5);

    // define the radius of the tool (sphere)
    double toolRadius = 0.01;

    // creates the tool
    tool = new cShapeSphere(toolRadius);
    world->addChild(tool);
    proxy = new cShapeSphere(toolRadius);
    world->addChild(proxy);
    proxy->setLocalPos(tool->getLocalPos());

    // set last positions
    xpbd_mesh->setVerticesLast();

    // get max vertex
    auto pos = xpbd_mesh->positions();
    vector<bool> indices(xpbd_mesh->numVerts());

    for (int i = 0; i < xpbd_mesh->numVerts() ; i++)
    {
        double height = pos(i,2);

        if (height > 0.300)
        {
            indices.at(i) = true;

        }
        else{
            indices.at(i) = false;
        }

    }

    xpbd_mesh->isFixed(indices);
    xpbd_mesh->createAABBCollisionDetector(toolRadius);


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
        //solve(xpbd_mesh, tool, Eigen::MatrixXd::Zero(xpbd_mesh->numVerts(),
        //3),0.01,20,5,true);

        timestep(xpbd_mesh, Eigen::MatrixXd::Zero(xpbd_mesh->numVerts(),3),0.005,50,5,true);


        updateGraphics();

        freqCounterGraphics.signal(1);

        // swap buffers
        glfwSwapBuffers(window);

        // process events
        glfwPollEvents();


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

    hapticDevice->close();

    // delete resources
    //delete hapticsThread;
    delete world;
    delete handler;

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

    // simulation in now running
    simulationRunning  = true;
    simulationFinished = false;

    // main haptic simulation loop
    while(simulationRunning) {
        /////////////////////////////////////////////////////////////////////
        // HAPTIC FORCE COMPUTATION
        /////////////////////////////////////////////////////////////////////

        // compute global reference frames for each object
        world->computeGlobalPositions(true);

        //cVector3d force(0,0,0);
        //hapticDevice->setForce(force);

        cVector3d tool_pos = tool->getLocalPos();
        cVector3d proxy_pos = proxy->getLocalPos();
        cVector3d n_ = cSub(tool->getLocalPos(), proxy->getLocalPos());
        double d_ = tool_pos.distance(proxy_pos);
        n_.normalize();
        double k = 1;

        //print force
        //std::cout << - n_.eigen().transpose()*d_*k << std::endl;

        // signal frequency counter
        freqCounterHaptics.signal(1);
    }

    // exit haptics thread
    simulationFinished = true;
}

void createTetrahedralMesh(cXPBDDeformableMesh* a_xpbdMesh)
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
            //a_chai3dMesh->newVertex(point);

            points.row(p) = Eigen::RowVector3d(output.pointlist[pi + 0],
                                               output.pointlist[pi + 1],
                                               output.pointlist[pi + 2]);

        }

        // sets the vertices of the mesh
        a_xpbdMesh->setVertices(points);
        a_xpbdMesh->buildAABBBoundaryBox(points);

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
    mass *= .00001;
    a_xpbdMesh->setMass(mass);
}

void timestep(            cXPBDDeformableMesh* model,
                          Eigen::MatrixX3d const& fext,
                          double timestep,
                          std::uint32_t iterations,
                          std::uint32_t substeps,
                          bool gravityEnabled)
{
    auto const num_iterations = iterations / substeps;
    double dt                 = timestep / static_cast<double>(substeps);
    auto const& constraints   = model->constraints();
    auto const J              = constraints.size();
    std::vector<double> lagrange_multipliers(J, 0.);

    // signal frequency counter
    auto fixed_ = model->fixed();
    proxy_stiffness = 10;

    for (auto s = 0u; s < substeps; ++s)
    {
        auto& v = model->velocity();
        auto& x = model->positions();

        auto const& m = model->mass();
        Eigen::MatrixX3d a = fext.array().colwise() / m.array();

        if (gravityEnabled)
        {
            Eigen::RowVector3d g(0,0,-9.81);
            a.rowwise() += g;
        }

        // explicit euler step
        auto vexplicit = v + dt * a;
        Eigen::MatrixXd p = x + dt * vexplicit;

        cVector3d goal_pos_temp;
        hapticDevice->getPosition(goal_pos_temp);

        auto goal_pos = goal_pos_temp.eigen();
        auto proxy_pos = proxy->getLocalPos().eigen();

        double t_;

        // Finds the collisions
        const auto& collisions = findCollisions(p,goal_pos,proxy_pos,*model,t_);
        if (collisions.collision == true)
        {
            computeCollisionConstraintsPassive(goal_pos, proxy_pos, p, x, model, t_, collisions);
            // std::cout << "Collision!" << std::endl;
        }
        else {
            proxy_pos = goal_pos;
        }

        tool->setLocalPos(goal_pos);
        proxy->setLocalPos(proxy_pos);

        // sequential gauss seidel type solve
        std::fill(lagrange_multipliers.begin(), lagrange_multipliers.end(), 0.0);

        for (auto n = 0u; n < num_iterations; ++n)
        {
            for (auto j = 0u; j < J; ++j)
            {
                auto const& constraint = constraints[j];
                constraint->project(p, v, x, m, lagrange_multipliers[j], dt);
            }
        }

        double lagrange_sum = 0;

        for (auto it = lagrange_multipliers.begin() ; it != lagrange_multipliers.end(); it++)
        {
            lagrange_sum += *it;
        }

        // update solution
        for (auto i = 0u; i < x.rows(); ++i)
        {
            if (fixed_.at(i) != true)
            {
                v.row(i) = (p.row(i) - x.row(i)) / dt;
                x.row(i) = p.row(i);
            }
        }

        model->updateChai3d(x);
    }
}

void computeCollisionConstraintsPassive(Eigen::Vector3d goal_ , Eigen::Vector3d& proxy_, Eigen::MatrixXd& p_, Eigen::MatrixXd& plast_,
                                        cXPBDDeformableMesh* a_mesh, double t_, const ColInfo& collisions)
{
    auto faces = a_mesh->faces();

    for (auto &col : collisions.faceCollisions)
    {
        auto face = faces.row(col);
        Eigen::Vector3d p1_ , p2_, p3_, p0_, n_;
        p1_ = p_.row(face(0));
        p2_ = p_.row(face(1));
        p3_ = p_.row(face(2));
        p0_ = (p1_ + p2_ + p3_) / 3 ;
        n_ = (p2_-p1_).cross(p3_-p1_).normalized();
        double A , B , C , D, dproxy_, dgoal_;
        A = n_(0); B = n_(1); C = n_(2);
        D = -A*p0_(0) - B*p0_(1) - C*p0_(2);
        dproxy_ = abs(A*proxy_(0) + B*proxy_(1) + C*proxy_(2) + D) / sqrt(pow(A,2) + pow(B,2) + pow(C,2));
        dgoal_ = abs(A*goal_(0) + B*goal_(1) + C*goal_(2) + D) / sqrt(pow(A,2) + pow(B,2) + pow(C,2));
        double breath_ = 0.001;
        Eigen::Vector3d proxyproj_ , goalproj_;

        if ((goal_ - proxy_).dot(n_) > 0)
        {
            // Object constraint
            // p_.row(face(0)) -= (d_ + breath_) * n_;
            // p_.row(face(1)) -= (d_ + breath_)* n_;
            // p_.row(face(2)) -= (d_ + breath_)* n_;

            // Proxy constraint
            proxyproj_ = proxy_ + dproxy_ * n_;
            goalproj_ = goal_ - dgoal_ * n_;

            proxy_ = goalproj_ - breath_ * n_;

        }
        else
        {
            // Object constraint
            // p_.row(face(0)) += (d_ + breath_) * n_;
            // p_.row(face(1)) += (d_ + breath_) * n_;
            // p_.row(face(2)) += (d_ + breath_) * n_;

            // Proxy constraint
            proxyproj_ = proxy_ - dproxy_ * n_;
            goalproj_ = goal_ + dgoal_ * n_;

            proxy_ = goalproj_ + breath_ * n_;
        }

        // p_.row(face(0)) = (p_.row(face(0)) - plast_.row(face(0))) * t_;
        // p_.row(face(1)) =  (p_.row(face(1)) - plast_.row(face(1))) * t_;
        // p_.row(face(2)) =  (p_.row(face(2)) - plast_.row(face(2))) * t_;

        //Eigen::Vector3d proxy_vel = (goal_ - proxy_) / dt;
    }
}

void computeCollisionConstraintsActive(Eigen::Vector3d goal_ , Eigen::Vector3d& proxy_, Eigen::MatrixXd& p_, Eigen::MatrixXd& plast_,
                                       cXPBDDeformableMesh* a_mesh, double t_, const ColInfo& collisions)
{

}

void testFriction(Eigen::Vector3d goal_ , Eigen::Vector3d& proxy_, Eigen::MatrixXd& p_, Eigen::MatrixXd& plast_,
                  cXPBDDeformableMesh* a_mesh, double t_, const ColInfo& collisions)
{

}
//--------------------