#ifndef SIMULATION_WIDGET
#define SIMULATION_WIDGET

// for Qt
#include "RGBABuffer.h"
#include "Vector2.h"
#include <QOpenGLWidget>
#include <vector>

// implementation of SPH based fluid simulation in 2D (Müller et al, SCA, 2003)

// a grid cell struct, contains references to its neighbours and a vector of particles
struct GridCell {
    // initialise the grid's neighbours to null, changed when generating a grid 
    GridCell() : cellParticles(0), neighbours(0) {}
    // a vector containing the ids of the particles that are in the cell. 
    // This is recomputed at each iteration, based on particle position (distance to origin)
    // we can determine the grid cell it belongs to with integer division
    std::vector<unsigned int> cellParticles;
    // pointers to the neighbours (top, left, bot, right) of the gridcell for searching particles
    std::vector<GridCell*> neighbours;
};

// struct container for a particle
struct Particle {
    Particle() : mass(0.5) {} // set mass to 0.5, vector default constructor is 0
    float mass;
    Vector2 position;
    Vector2 velocity;
    GridCell* cell; // some unique id referring to a gridcell
};


class SimulationWidget : public QOpenGLWidget {
    Q_OBJECT

    public:
    // constructor
    SimulationWidget(QWidget* parent);
    // destructor
    ~SimulationWidget();

    public slots:
    // slots for updating the simulation parameters
    void updateParticleCount(int newCount);
    void updateParticleRadius(int newRadius);
    void updateBlobHeight(int newHeight);
    void updateGravityConstant(int newGravity);
    void updateRestDensity(int newDensity);
    void updateGasConstant(int newConstant);
    void updateViscosityConstant(int newConstant);
    void updateSurfaceTensionConstant(int newConstant);
    void updateCoreRadius(int newRadius);
    void usePoly6();
    void useSpiky();
    void addViscosity(int state);
    void addSurfaceTension(int state);

    // called whenever the user presses rest button
    void resetSimulation();

    // intergrate and paint
    void stepSimulation();

    protected:
    //
    // Qt opengl functions
    //

	// called when OpenGL context is set up
	void initializeGL();
	// called every time the widget is resized
	void resizeGL(int w, int h);
	// called every time the widget needs painting
	void paintGL();

    private:
    //
    // Scene generation
    //

    // generates the grid
    void GenerateGrid();
    // generate the particles
    void GenerateParticles();
    // assigns particles to their corresponding grid cell
    void AssignCellParticles();

    //
    // Forces computation and integration
    //

    // compute the density of a particle
    float ComputeDensity(const unsigned int& particle);
    // compute the pressure force applied on a particle
    Vector2 ComputePressurePoly6(const unsigned int& particle, const float& p_i);
    Vector2 ComputePressureSpiky(const unsigned int& particle, const float& p_i);
    // compute the viscosity force applied on a particle
    Vector2 ComputeViscosity(const unsigned int& particle);
    // compute the surface tension force applied on a particle
    Vector2 ComputeSurfaceTension(const unsigned int& particle);
    // return a vector modified by collision (if any)
    void ComputeCollision(const unsigned int& particle);
    // use leapfrog integration
    void LeapFrogIntegrate();

    //
    // Drawing routines
    //

    // draws the grid cells on the framebuffer
    void DrawGridCells();
    // draws a small circle for each particle
    void DrawParticles();
    // colours pixels to show the tank
    void DrawTank();

    private:
    // pointer to an image buffer where the simulation is drawn
    RGBABuffer frameBuffer;

    // a vector containing all the particles (of size particleCount)
    std::vector<Particle> particles;

    // a 2D vector containing all the grid cells
    std::vector<std::vector<GridCell>> grid;

    // some scene settings
    const float sceneSize;
    float aspectRatio;
    int particleRadius;
    // size (length and width) of the rectangular grid
    int gridSize;

    // variables that are set in the UI:
    // number of particles in the scene
    int particleCount;
    // gravity constant
    float gravityConstant;
    // rest density
    float rhoRest;
    // gas constant k
    float gasK;
    // viscosity mu
    float viscosityMu;
    // surface tension sigma
    float tensionSigma;
    // kernel core radius h
    float coreRadius;
    // height of the water blob
    float blobHeight;
    // aribtrary time step for integration 
    const float deltaT;
    bool usingPoly6;
    bool addingViscosity;
    bool addingSurfaceTension;
};


#endif