// the window class declaration
#include "Window.h"

// constructor for control button
CtrlButton::CtrlButton(QString label, QWidget *parent) 
    : QPushButton(label, parent) {
    setFixedSize(70, 30);
}

// constructor
Window::Window(QWidget* parent) : QWidget(parent) {

    //
    // SET UP THE WIDGETS
    //

    // give the window a good name
    setWindowTitle("Simple Fluid Sim");
    // create the simulation widget and set the window as its parent
    simulationWidget = new SimulationWidget(this);
    simulationWidget->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);

    // create the timer widget and connect it
    timer = new QTimer(this);
    timer->setInterval(16); // 1 / 60 seconds
    // upon timeout, update the scene
    QObject::connect(timer, SIGNAL(timeout()), simulationWidget, SLOT(stepSimulation())); 

    //
    // SETUP THE PLAYBACK 
    // 

    // init simulation play controls
    playbackLayout = new QHBoxLayout;
    // init simulation buttons
    play = new CtrlButton("&play", this);
    stop = new CtrlButton("&stop", this);
    reset = new CtrlButton("&reset", this);
    // connect to timer
    QObject::connect(play, SIGNAL(pressed()), timer, SLOT(start()));
    QObject::connect(stop, SIGNAL(pressed()), timer, SLOT(stop()));
    QObject::connect(reset, SIGNAL(pressed()), timer, SLOT(stop()));
    QObject::connect(reset, SIGNAL(pressed()), 
        simulationWidget, SLOT(resetSimulation()));
    // add to the control layout
    playbackLayout->addWidget(play);
    playbackLayout->addWidget(stop);
    playbackLayout->addWidget(reset);

    //
    // SETUP THE PARAMETERS 
    //

    QSlider* particleRadiusSlider;
    QLabel* particleRadiusSliderLabel;

    // init the layout and group box
    parametersGroupBox = new QGroupBox("Parameters");
    parametersLayout = new QVBoxLayout;
    // create the widgets (mostly sliders)
    particleCountSliderLabel = new QLabel("Number of particles (100 to 10000)", this);
    particleCountSlider = new QSlider(Qt::Horizontal, this);
    particleRadiusSliderLabel = new QLabel("Radius of a particle (5 to 20)", this);
    particleRadiusSlider = new QSlider(Qt::Horizontal, this);
    heightSliderLabel  = new QLabel("Initial height (2 to 6)", this);
    heightSlider = new QSlider(Qt::Horizontal, this);
    gravitySliderLabel  = new QLabel("Gravity (0 to 10)", this);
    gravitySlider = new QSlider(Qt::Horizontal, this);
    restDensitySliderLabel  = new QLabel("Rest density (0 to 1)", this);
    restDensitySlider = new QSlider(Qt::Horizontal, this);
    gasKSliderLabel  = new QLabel("Gas constant (0 to 1)", this);
    gasKSlider = new QSlider(Qt::Horizontal, this);
    viscosityMuSliderLabel  = new QLabel("Viscosity constant (0 to 5)", this);
    viscosityMuSlider = new QSlider(Qt::Horizontal, this);
    tensionSigmaSliderLabel  = new QLabel("Surface tension constant (0 to 1)", this);
    tensionSigmaSlider = new QSlider(Qt::Horizontal, this);
    coreRadiusSliderLabel  = new QLabel("Core radius (0 to 1)", this);
    coreRadiusSlider = new QSlider(Qt::Horizontal, this);
    pressureButtonGroup = new QButtonGroup;
    poly6PressureBox = new QCheckBox("Poly6 pressure");
    spikyPressureBox = new QCheckBox("Spiky pressure");
    viscosityBox = new QCheckBox("Viscosity");
    surfaceTensionBox = new QCheckBox("Surface tension");
    // connect the widgets
    QObject::connect(particleCountSlider, SIGNAL(valueChanged(int)), 
        simulationWidget, SLOT(updateParticleCount(int)));
    QObject::connect(particleRadiusSlider, SIGNAL(valueChanged(int)), 
        simulationWidget, SLOT(updateParticleRadius(int)));
    QObject::connect(heightSlider, SIGNAL(valueChanged(int)), 
        simulationWidget, SLOT(updateBlobHeight(int)));
    QObject::connect(gravitySlider, SIGNAL(valueChanged(int)), 
        simulationWidget, SLOT(updateGravityConstant(int)));
    QObject::connect(restDensitySlider, SIGNAL(valueChanged(int)), 
        simulationWidget, SLOT(updateRestDensity(int)));
    QObject::connect(gasKSlider, SIGNAL(valueChanged(int)), 
        simulationWidget, SLOT(updateGasConstant(int)));
    QObject::connect(viscosityMuSlider, SIGNAL(valueChanged(int)), 
        simulationWidget, SLOT(updateGasConstant(int)));
    QObject::connect(tensionSigmaSlider, SIGNAL(valueChanged(int)), 
        simulationWidget, SLOT(updateSurfaceTensionConstant(int)));
    QObject::connect(coreRadiusSlider, SIGNAL(valueChanged(int)), 
        simulationWidget, SLOT(updateCoreRadius(int)));
    QObject::connect(poly6PressureBox, SIGNAL(pressed()), 
        simulationWidget, SLOT(usePoly6()));
    QObject::connect(spikyPressureBox, SIGNAL(pressed()), 
        simulationWidget, SLOT(useSpiky()));
    QObject::connect(viscosityBox, SIGNAL(stateChanged(int)), 
        simulationWidget, SLOT(addViscosity(int)));
    QObject::connect(surfaceTensionBox, SIGNAL(stateChanged(int)), 
        simulationWidget, SLOT(addSurfaceTension(int)));
    // set slider parameters
    particleCountSlider->setRange(1, 10); // slider scaled by 10
    particleCountSlider->setValue(1);
    particleRadiusSlider->setRange(5, 20); // not scaled
    particleRadiusSlider->setValue(10);
    heightSlider->setRange(200, 600); // slider scaled by 1/100
    heightSlider->setValue(400);
    gravitySlider->setRange(0, 100); // slider scaled by 1/10
    gravitySlider->setValue(98);
    restDensitySlider->setRange(0, 100); // slider scaled by 1/100
    restDensitySlider->setValue(0);
    gasKSlider->setRange(1, 100); // slider scaled by 1/100
    gasKSlider->setValue(50);
    viscosityMuSlider->setRange(1, 50); // slider scaled by 1/10
    viscosityMuSlider->setValue(25);
    tensionSigmaSlider->setRange(1, 100); // slider scaled by 1/100
    tensionSigmaSlider->setValue(100); 
    coreRadiusSlider->setRange(1, 100); // slider  scaled by 1/100
    coreRadiusSlider->setValue(100);
    pressureButtonGroup->addButton(poly6PressureBox);
    pressureButtonGroup->addButton(spikyPressureBox);
    pressureButtonGroup->setExclusive(true);
    poly6PressureBox->setCheckState(Qt::Checked);
    // set maximum group size
    parametersGroupBox->setMaximumWidth(300);
    // add them to the layout
    parametersLayout->addWidget(particleCountSliderLabel);
    parametersLayout->addWidget(particleCountSlider);
    parametersLayout->addWidget(particleRadiusSliderLabel);
    parametersLayout->addWidget(particleRadiusSlider);
    parametersLayout->addWidget(heightSliderLabel);
    parametersLayout->addWidget(heightSlider);
    parametersLayout->addWidget(gravitySliderLabel);
    parametersLayout->addWidget(gravitySlider);
    parametersLayout->addWidget(restDensitySliderLabel);
    parametersLayout->addWidget(restDensitySlider);
    parametersLayout->addWidget(gasKSliderLabel);
    parametersLayout->addWidget(gasKSlider);
    parametersLayout->addWidget(viscosityMuSliderLabel);
    parametersLayout->addWidget(viscosityMuSlider);
    parametersLayout->addWidget(coreRadiusSliderLabel);
    parametersLayout->addWidget(coreRadiusSlider);
    parametersLayout->addWidget(tensionSigmaSliderLabel);
    parametersLayout->addWidget(tensionSigmaSlider);
    parametersLayout->addWidget(poly6PressureBox);
    parametersLayout->addWidget(spikyPressureBox);
    parametersLayout->addWidget(viscosityBox);
    parametersLayout->addWidget(surfaceTensionBox);
    // set the group box's layout
    parametersGroupBox->setLayout(parametersLayout);

    //
    // SET UP THE WINDOW'S LAYOUT
    // 

    // init the layouts
    mainLayout = new QHBoxLayout(this);
    simulationLayout = new QVBoxLayout;
    // add the menu bar to the layout
    // add the sub layouts
    mainLayout->addLayout(simulationLayout);
    mainLayout->addWidget(parametersGroupBox);
    // add the simulation widget to the simulation layout
    simulationLayout->addWidget(simulationWidget);
    // add the simulation controls to the layout
    simulationLayout->addLayout(playbackLayout);

    // reset the simulation to initiaite scene setup
    simulationWidget->resetSimulation();
}

// destructor
Window::~Window() {
    // do something maybe?
}
