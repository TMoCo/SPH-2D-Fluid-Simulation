#ifndef WINDOW_H
#define WINDOW_H

// standard Qt widgets
#include <QtWidgets>

// the simulation widget
#include "SimulationWidget.h"

// convenience class for animation control buttons, has a set width and height
class CtrlButton : public QPushButton {
    Q_OBJECT
    public:
    //  constructor, sets the geometry of button 
    CtrlButton(QString label, QWidget *parent);
};

// Window class that behaves as a root widget for the application
class Window : public QWidget {
    Q_OBJECT

    public:
    // constructor
    Window(QWidget* parent);
    // destructor
    ~Window();

    public slots:
    
    // signals for file loading and saving
    signals:


    private:
    //
    // MAIN WIDGET AND LAYOUTS
    //

    // window main layout
    QHBoxLayout* mainLayout;
    // layout for simulation widgets (openGL widget and player)
    QVBoxLayout* simulationLayout;

    // widget where the simulation is rendered
    SimulationWidget* simulationWidget;

    //
    // SIMULATION PLAYBACK
    //

    // timer for updating the scene
    QTimer* timer;

    // layout for the simulation player
    QHBoxLayout* playbackLayout;
    // widgets for controlling the simulation
    CtrlButton* play;
    CtrlButton* stop;
    CtrlButton* reset;

    //
    // SIMULATION PARAMETERS
    //

    // parameters group box 
    QGroupBox* parametersGroupBox;
    QVBoxLayout* parametersLayout;

    // sliders
    // sets the number of particles in the scene
    QSlider* particleCountSlider;
    QLabel* particleCountSliderLabel;
    // sets the radius of a particle
    QSlider* particleRadiusSlider;
    QLabel* particleRadiusSliderLabel;
    // modifies the initial height of the fluid blob
    QSlider* heightSlider;
    QLabel* heightSliderLabel;
    // modifies the intensity of the gravity
    QSlider* gravitySlider;
    QLabel* gravitySliderLabel;
    // modifies the rest density
    QSlider* restDensitySlider;
    QLabel* restDensitySliderLabel;
    // gas constant k
    QSlider* gasKSlider;
    QLabel* gasKSliderLabel;
    // viscosity constant mu
    QSlider* viscosityMuSlider;
    QLabel* viscosityMuSliderLabel;
    // tension constant sigma
    QSlider* tensionSigmaSlider;
    QLabel* tensionSigmaSliderLabel;
    // number of cells slider
    QSlider* coreRadiusSlider;
    QLabel* coreRadiusSliderLabel;

    // check boxes
    // checkboxes for switching options on or off
    QCheckBox* poly6PressureBox; 
    QCheckBox* spikyPressureBox;
    QCheckBox* viscosityBox;
    QCheckBox* surfaceTensionBox;
    // grouping of the pressure check boxes to make them mutually exclusive
    QButtonGroup* pressureButtonGroup;
};

#endif