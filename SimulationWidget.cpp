#include <math.h>

#include "Kernel.h"
#include "SimulationWidget.h"

// constructor
SimulationWidget::SimulationWidget(QWidget* parent) : QOpenGLWidget(parent), 
    // scene of size 10 (m) and delta t of 16 ms. Other variables will be set after widget construction
    sceneSize(10.0), deltaT(0.016), usingPoly6(true), addingViscosity(false), 
    addingSurfaceTension(false) {}

// destructor
SimulationWidget::~SimulationWidget() {}

//
// OPENGL METHODS
//

void SimulationWidget::initializeGL() {}

// called every time the widget is resized
void SimulationWidget::resizeGL(int w, int h) {
    // resize the framebuffer
    frameBuffer.Resize(w, h);
    aspectRatio = h / (float)w;
}
	
// called every time the widget needs painting
void SimulationWidget::paintGL() {
    // set background colour
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    // set the buffer to a light blue as a background colour
    for (long r = 0; r < frameBuffer.height; r++)
        for (long c = 0; c < frameBuffer.width; c++)
            frameBuffer[r][c] = RGBAValue<float>(255, 225, 255, 255);

    // put the scene in the framebuffer
    DrawTank();
    DrawParticles();

    // and display the contents of the framebuffer
    glDrawPixels(frameBuffer.width, frameBuffer.height, 
        GL_RGBA, GL_UNSIGNED_BYTE, frameBuffer.block);
}

//
// SCENE GENERATION
//

void SimulationWidget::GenerateGrid() {
    // start by erasing the previous grid
    grid.resize(0);
    // compute the grid dimension using the current core radius and scene dimension
    gridSize = sceneSize / coreRadius;
    // now resize the grid to accomodate the desired number of grid elements
    grid.resize(gridSize);
    for (size_t cell = 0; cell < grid.size(); cell++)
        grid[cell].resize(gridSize);
    // now loop over each grid element and assign its neighbours depending on the index
    for (int row = 0; row < gridSize; row++)
        for (int col = 0; col < gridSize; col++) {
            // depending on the id, update the neighbour pointer
            if (row != 0) // has a left neighbour
                grid[row][col].neighbours.push_back(&grid[row-1][col]);
            if (row != gridSize-1) // has a right neighbour
                grid[row][col].neighbours.push_back(&grid[row+1][col]);
            if (col != 0) // has a bottom neighbour
                grid[row][col].neighbours.push_back(&grid[row][col-1]);
            if (col != gridSize-1) // has a top neighbour
                grid[row][col].neighbours.push_back(&grid[row][col+1]);
        }
}

// generate the particles in the water blob
void SimulationWidget::GenerateParticles() {
    // start by erasing previous particles
    particles.resize(0);
    // resize to accomodate the right number of particles
    particles.resize(particleCount*particleCount);
    // initial position in space at bottom left of the blob
    float startX = 2;
    float startY = blobHeight;
    //unsigned int id = 0;
    // loop over the particles to generate
    for (int i = 0; i < particleCount; i++)
        for (int j = 0; j < particleCount; j++)
            particles[i*particleCount+j].position =
                Vector2(startX + i * 2 / (float)particleCount, startY + j * 3 / (float)particleCount);
}

// assigns particles to their corresponding grid cell (in world space!!!)
void SimulationWidget::AssignCellParticles() {
    // reset the grid
    for (int row = 0; row < gridSize; row++)
        for (int col = 0; col < gridSize; col++)
            grid[row][col].cellParticles.resize(0);
    // loop over the scene particles and determine which grid cell they belong to
    // using modular arithmetic
    int row, col;
    for (unsigned int p = 0; p < particles.size(); p++) {
        // get the particle's cell by dividing by the core radius, cast to int
        row = particles[p].position.y / coreRadius;
        col = particles[p].position.x / coreRadius;
        // determine which cell the pixel is in
        grid[row][col].cellParticles.push_back(p);
        particles[p].cell = &grid[row][col];
    }
}

//
// FORCE COMPUTATION AND INTEGRATION
//

// compute the density of a given particle
float SimulationWidget::ComputeDensity(const unsigned int& particle) {
    // initialise the sum at 0
    Vector2 distance; 
    float density = 0;
    // loop around the particles in the grid cell the particle belongs to
    GridCell* cell = particles[particle].cell;
    for (size_t p = 0; p < cell->cellParticles.size(); p++) {
        // skip the current particle
        //if (cell->cellParticles[p] == particle)
        //    continue;
        // get the vector from particle to other
        distance = particles[particle].position - 
            particles[cell->cellParticles[p]].position;
        // sum up 
        density += particles[cell->cellParticles[p]].mass * Kernel::Poly6(distance, coreRadius);
    }
    // do the same for the neighbouring grid cells
    for (size_t neighbour = 0; neighbour < cell->neighbours.size(); neighbour++ ) {
        // loop over the neighbour's particles, no need to check for current particle
        for (size_t p = 0; p < cell->neighbours[neighbour]->cellParticles.size(); p++) {
            // get the vector from particle to other
            distance = particles[particle].position - 
                particles[cell->neighbours[neighbour]->cellParticles[p]].position;
            // sum up 
            density += particles[cell->neighbours[neighbour]->cellParticles[p]].mass * 
                Kernel::Poly6(distance, coreRadius);
        }
    }
    return density;
}

// computes and returns the pressure force for a single particle
Vector2 SimulationWidget::ComputePressurePoly6(const unsigned int& particle, const float& p_i) {
    // initialise pressure force at 0
    Vector2 distance, forceDirection;
    float pressure = 0, rho_j;
    // loop around the particles in the grid cell the particle belongs to
    GridCell* cell = particles[particle].cell;
    for (size_t p = 0; p < cell->cellParticles.size(); p++) {
        // skip the current particle
        //if (cell->cellParticles[p] == particle)
        //    continue;
        rho_j = ComputeDensity(cell->cellParticles[p]);
        // skip bad density
        if(rho_j <= EPSILON)
            continue;
        // get the vector from particle to other
        distance = particles[particle].position - particles[cell->cellParticles[p]].position;
        // sum up 
        forceDirection += distance;
        pressure += particles[cell->cellParticles[p]].mass * 
            ((p_i + gasK * rho_j) / (2 * (rho_j- rhoRest))) * Kernel::Poly6Gradient(distance, coreRadius);
    }
    // do the same for the neighbouring grid cells
    for (size_t neighbour = 0; neighbour < cell->neighbours.size(); neighbour++ )
        // loop over the neighbour's particles, no need to check for current particle
        for (size_t p = 0; p < cell->neighbours[neighbour]->cellParticles.size(); p++) {
            rho_j = ComputeDensity(cell->neighbours[neighbour]->cellParticles[p]);
            // skip bad density
            if(rho_j <= EPSILON)
                continue;
            // get the vector from particle to other
            distance = particles[particle].position - 
                particles[cell->neighbours[neighbour]->cellParticles[p]].position;
            // sum up 
            forceDirection += distance;
            pressure += particles[cell->neighbours[neighbour]->cellParticles[p]].mass *
                ((p_i + gasK * rho_j) / (2 * (rho_j- rhoRest))) * Kernel::Poly6Gradient(distance, coreRadius);
        }
    return -pressure * forceDirection.unit();//
}

// computes and returns the pressure force for a single particle
Vector2 SimulationWidget::ComputePressureSpiky(const unsigned int& particle, const float& p_i) {
    // initialise pressure force at 0
    Vector2 distance, forceDirection;
    float pressure = 0, rho_j;
    
    // loop around the particles in the grid cell the particle belongs to
    GridCell* cell = particles[particle].cell;
    for (size_t p = 0; p < cell->cellParticles.size(); p++) {
        // skip the current particle
        //if (cell->cellParticles[p] == particle)
        //    continue;
        rho_j = ComputeDensity(cell->cellParticles[p]);
        // skip bad density
        if(rho_j <= EPSILON)
            continue;
        // get the vector from particle to other
        distance = particles[particle].position - particles[cell->cellParticles[p]].position;
        // sum up 
        forceDirection += distance;
        pressure += particles[cell->cellParticles[p]].mass *
            ((p_i + gasK * rho_j) / (2 * (rho_j- rhoRest))) * Kernel::SpikyGradient(distance, coreRadius);
    }
    // do the same for the neighbouring grid cells
    for (size_t neighbour = 0; neighbour < cell->neighbours.size(); neighbour++ )
        // loop over the neighbour's particles, no need to check for current particle
        for (size_t p = 0; p < cell->neighbours[neighbour]->cellParticles.size(); p++) {
            rho_j = ComputeDensity(cell->neighbours[neighbour]->cellParticles[p]);
            // skip bad density
            if(rho_j <= EPSILON)
                continue;
            // get the vector from particle to other
            distance = particles[particle].position - 
                particles[cell->neighbours[neighbour]->cellParticles[p]].position;
            // sum up 
            forceDirection += distance;
            pressure += particles[cell->neighbours[neighbour]->cellParticles[p]].mass *
                ((p_i + gasK * rho_j) / (2 * (rho_j- rhoRest))) * Kernel::SpikyGradient(distance, coreRadius);
        }
    return -pressure * forceDirection.unit();
}

// computes and returns the viscosity force for a single particle
Vector2 SimulationWidget::ComputeViscosity(const unsigned int& particle) {
    Vector2 viscosity, particleVelocity = particles[particle].velocity;
    float rho_j;

    GridCell* cell = particles[particle].cell;
    for (size_t p = 0; p < cell->cellParticles.size(); p++) {
        rho_j = ComputeDensity(cell->cellParticles[p]);
        // skip bad density
        if(rho_j <= EPSILON)
            continue;
        // sum up 
        viscosity += particles[p].mass * ((particles[cell->cellParticles[p]].velocity - particleVelocity) / rho_j) * 
            Kernel::ViscosityLaplacian(particles[particle].position - particles[cell->cellParticles[p]].position, coreRadius);
    }
    // do the same for the neighbouring grid cells
    for (size_t neighbour = 0; neighbour < cell->neighbours.size(); neighbour++ ) {
        // loop over the neighbour's particles, no need to check for current particle
        for (size_t p = 0; p < cell->neighbours[neighbour]->cellParticles.size(); p++) {
            rho_j = ComputeDensity(cell->neighbours[neighbour]->cellParticles[p]);
            // skip bad density
            if(rho_j <= EPSILON)
                continue;
            // sum up 
            viscosity += particles[p].mass * ((particles[cell->neighbours[neighbour]->cellParticles[p]].velocity - 
                particleVelocity) / rho_j) * Kernel::ViscosityLaplacian(particles[particle].position - 
                particles[cell->neighbours[neighbour]->cellParticles[p]].position, coreRadius);
        }
    }
    return viscosityMu * viscosity;
}

// computes and returns the surface tension for a single particle
Vector2 SimulationWidget::ComputeSurfaceTension(const unsigned int& particle) {
    Vector2 distance, forceDirection, normal;
    float colourLaplacian = 0, rho_j;

    GridCell* cell = particles[particle].cell;
    for (size_t p = 0; p < cell->cellParticles.size(); p++) {
        // skip the current particle
        if (cell->cellParticles[p] == particle)
            continue;
        rho_j = ComputeDensity(cell->cellParticles[p]);
        // skip bad density
        if (rho_j <= EPSILON)
            continue;
        // get the vector from particle to other
        distance = particles[particle].position - particles[cell->cellParticles[p]].position;
        // sum up
        forceDirection += distance;
        // normal 
        normal += particles[cell->cellParticles[p]].mass * (1.0 / rho_j) *
            Kernel::Poly6Gradient(distance, coreRadius) * distance;
        // kappa
        colourLaplacian += particles[cell->cellParticles[p]].mass * (1.0 / rho_j) *
            Kernel::Poly6Laplacian(distance, coreRadius);
    }
    // do the same for neighbouring grid cells
    for (size_t neighbour = 0; neighbour < cell->neighbours.size(); neighbour++)
            for (size_t p = 0; p < cell->neighbours[neighbour]->cellParticles.size(); p++) {
                rho_j = ComputeDensity(cell->neighbours[neighbour]->cellParticles[p]);
                // skip bad density
                if (rho_j <= EPSILON)
                    continue;
                // get the vector from particle to other
                distance = particles[particle].position - 
                    particles[cell->neighbours[neighbour]->cellParticles[p]].position;
                // sum up
                forceDirection += distance;
                // normal 
                normal += particles[cell->neighbours[neighbour]->cellParticles[p]].mass *
                    (1.0 / rho_j) * Kernel::Poly6Gradient(distance, coreRadius) * distance;
                // kappa
                colourLaplacian += particles[cell->neighbours[neighbour]->cellParticles[p]].mass *
                    (1.0 / rho_j) * Kernel::Poly6Laplacian(distance, coreRadius);
        }
    // normal
    float normalLength = normal.length();
    if (normalLength < EPSILON) 
        return Vector2();
    else 
        return -tensionSigma * (colourLaplacian / normalLength) * normal.unit();
}
// use leapfrog integration to update the scene
void SimulationWidget::LeapFrogIntegrate() {
    Vector2 v_ihalf, forces;
    float rho_i;
    // leap frog for each particle 
    // loop over the particles to assign the new positions
    for (size_t p = 0; p < particles.size(); p++) {
        // velocity + acceleration * 1/2dt
        v_ihalf = particles[p].velocity + Vector2(0, -gravityConstant) * (deltaT / 2.0);
        // start by adding external forces
        // get the density
        rho_i = ComputeDensity(p);
        // only add the forces to the velocity if the density is valid
        if (rho_i > EPSILON) {
            forces.x = forces.y = 0;
            // depending on UI, compute forces
            if (usingPoly6)
                forces += ComputePressurePoly6(p, gasK * (rho_i - rhoRest)) / (rho_i) ;
            else
                forces += ComputePressureSpiky(p, gasK * (rho_i - rhoRest)) / (rho_i) ;
            if (addingViscosity)
                forces += ComputeViscosity(p) / (rho_i) ;
            if (addingSurfaceTension)
                forces += ComputeSurfaceTension(p) / (rho_i) ;
            // now compute velocity at half time step
            v_ihalf += forces * (deltaT / 2.0);
        }
        //else ;
        // p_n+1 = p_n + v_ihalf * dt
        particles[p].position = particles[p].position + v_ihalf * deltaT;
        // save half velocity for now
        particles[p].velocity = v_ihalf;
    }
    // sum half velocity with new half acceleration
    for (size_t p = 0; p < particles.size(); p++) {
        Vector2 a_ihalf = Vector2(0, -gravityConstant);
        float rho_i = ComputeDensity(p);
        // skip bad density
        // other wise add half velocity
        if (rho_i > EPSILON) {
            forces.x = forces.y = 0; 
            // depending on UI, compute forces
            if (usingPoly6)
                forces += ComputePressurePoly6(p, gasK * (rho_i - rhoRest)) / (rho_i) ;
            else
                forces += ComputePressureSpiky(p, gasK * (rho_i - rhoRest)) / (rho_i) ;
            if (addingViscosity)
                forces += ComputeViscosity(p) / (rho_i) ;
            if (addingSurfaceTension)
                forces += ComputeSurfaceTension(p) / (rho_i) ;
            // now compute velocity at half time step
            a_ihalf += forces;
        }
        //else 
        particles[p].velocity += a_ihalf * (deltaT / 2.0);
    }
    
    // compute collisions with the tank (and add a ceiling for convenience)
    float dampen = 0.5f;
    for (size_t p = 0; p < particles.size(); p++) {
        if (particles[p].position.y < 1) {
            particles[p].velocity.y = -particles[p].velocity.y * dampen;
            particles[p].position.y = 1;
        }
        if (particles[p].position.y > 9.5) {
            particles[p].velocity.y = -particles[p].velocity.y * dampen;
            particles[p].position.y = 9.5;
        }
        if (particles[p].position.x < 1) {
            particles[p].velocity.x = -particles[p].velocity.x * dampen;
            particles[p].position.x = 1;
        }
        if (particles[p].position.x > 9) {
            particles[p].velocity.x = -particles[p].velocity.x * dampen;
            particles[p].position.x = 9;
        }
    }
}

//
// DRAWING ROUTINES
//

// Debug routine that colours frame buffer pixels depending on grid cells 
void SimulationWidget::DrawGridCells() {
    // get the number of pixels per grid
    float cellHeight, cellWidth;
    if (aspectRatio < 1.0) {
        cellHeight = frameBuffer.height / sceneSize * coreRadius;
        cellWidth = frameBuffer.width / sceneSize * coreRadius * aspectRatio;
    }
    else {
        cellHeight = frameBuffer.height / sceneSize * coreRadius * aspectRatio;
        cellWidth = frameBuffer.width / sceneSize * coreRadius;
    }
    // colour for a cell that contins a particle
    RGBAValue<float> cellColour(0, 255, 0, 255);
    // loop over the grid
    for (int row = 0; row < gridSize; row++)
        for (int col = 0; col < gridSize; col++) {
            // for now, just paint non empty cell same green
            if (grid[row][col].cellParticles.empty())
                continue;
            // based on the row and col of each cell, colour the pixels. Use floats
            // to make sure we reach the whole frame buffer then cast to int
            if ((row * cellHeight > frameBuffer.height) || (col * cellWidth > frameBuffer.width))
                continue;
            for (float h = 0; h < cellHeight; h++) // row
                for (float w = 0 ; w < cellWidth; w++) { // col 
                    // set the colour here... we will use the grid cell's density
                    // to determine the cell pixel's colour
                    frameBuffer[(int)(row * cellHeight + h)][(int)(col * cellWidth + w)] = cellColour;
                }

        }
}

// Routine to draw individual particles
void SimulationWidget::DrawParticles() {
    // loop over the particles and draw at pixel location
    int particleRow, particleCol; // pixel representing particle position
    int blobRow, blobCol;
    RGBAValue<float> colour(255, 0, 0, 255);
    for (size_t p = 0; p < particles.size(); p++) {
        if (aspectRatio < 1.0) {
            particleRow = particles[p].position.y * frameBuffer.height / sceneSize;
            particleCol = particles[p].position.x * frameBuffer.width / sceneSize  * aspectRatio;;
        }
        else {
            particleCol = (particles[p].position.x * frameBuffer.width) / sceneSize;
            particleRow = particles[p].position.y * frameBuffer.height / sceneSize * aspectRatio;
        }
        // now draw a sphere around that pixel using a bounding box
        for (int i = 0; i < 2 * particleRadius; i++) {
            blobRow = particleRow - particleRadius + i;
            for (int j = 0; j < 2 * particleRadius; j++) {
                blobCol = particleCol - particleRadius + j;
                // vector from blob pixel to particle pixel
                Vector2 toParticle(particleCol - blobCol, particleRow - blobRow);
                if (toParticle.length() > particleRadius)
                    continue;
                frameBuffer[blobRow][blobCol] = colour;
            }
        }

    }
}

// determine the area to be coloured to represent a tank
void SimulationWidget::DrawTank() {
    // the tank will be 5 units wide and five units high with the bottom at 1. 
    RGBAValue<float> tankColour(51, 102, 153, 255);
    // the tank boundaries will be about 20 pixels thick
    float thick = 15;
    float tankHeight, tankWidth, startH, startW;
    if (aspectRatio < 1.0) {
        tankHeight = (5 * frameBuffer.height) / sceneSize;
        tankWidth = (8 * frameBuffer.width) / sceneSize * aspectRatio;
        startH = frameBuffer.height / sceneSize;
        startW = frameBuffer.width / sceneSize * aspectRatio;
    }
    else {
        tankHeight = (5 * frameBuffer.height) / sceneSize * aspectRatio;
        tankWidth = (8 * frameBuffer.width) / sceneSize;
        startH = frameBuffer.height / sceneSize * aspectRatio;
        startW = frameBuffer.width / sceneSize;
    }
    // left side
    // colour the left wall, starting at the bottom
    for (float h = 0; h < tankHeight; h++) // row
        for (float w = 0 ; w < thick + particleRadius; w++) // col 
            frameBuffer[(int)(startH + h)][(int)(startW - w - particleRadius)] = tankColour;
    // colour the bottom of the tank, starting at the left
    for (float h = 0; h < thick + particleRadius; h++) // row
        for (float w = 0 ; w < tankWidth; w++) // col 
            frameBuffer[(int)(startH - h - particleRadius)][(int)(startW + w)] = tankColour;
    // colour the right wall, starting at the bottom
    if (aspectRatio < 1.0) 
        startW = (9 * frameBuffer.width) / sceneSize * aspectRatio;
    else
        startW = (9 * frameBuffer.width) / sceneSize;    
    for (float h = 0; h < tankHeight; h++) // row
        for (float w = 0 ; w < thick + particleRadius; w++) // col 
            frameBuffer[(int)(startH + h)][(int)(startW + w + particleRadius)] = tankColour;
}

//
// SLOTS
//

// slots to update various scene variables
void SimulationWidget::updateParticleCount(int newCount) {
    // scaled by 10
    particleCount = newCount * 10;
}

void SimulationWidget::updateParticleRadius(int newRadius) {
    particleRadius = newRadius;
}

void SimulationWidget::updateBlobHeight(int newHeight) {
    // scaled by 1/100
    blobHeight = newHeight / 100.0;
}

void SimulationWidget::updateGravityConstant(int newGravity) {
    // scaled by 1/10
    gravityConstant = newGravity / 10.0;
}

void SimulationWidget::updateRestDensity(int newDensity) {
    // scaled by 100
    rhoRest = newDensity / 100.0;
}

void SimulationWidget::updateGasConstant(int newConstant) {
    // scaled by 1/100
    gasK = newConstant / 100.0;
}

void SimulationWidget::updateViscosityConstant(int newConstant) {
    // scaled by 1/100
    viscosityMu = newConstant / 100.0;
}

void SimulationWidget::updateSurfaceTensionConstant(int newConstant) {
    // scaled by 1/100
    tensionSigma = newConstant / 10.0;
}

void SimulationWidget::updateCoreRadius(int newRadius) {
    // scaled by 1/100
    coreRadius = newRadius / 100.0;
}

void SimulationWidget::usePoly6() {
    usingPoly6 = true;
}

void SimulationWidget::useSpiky() {
    usingPoly6 = false;
}

void SimulationWidget::addViscosity(int state) {
    addingViscosity = state;
}

void SimulationWidget::addSurfaceTension(int state) {
    addingSurfaceTension = state;
    std::cout << addingSurfaceTension << std::endl;
}


// reset the simulation
void SimulationWidget::resetSimulation() {
    // generate a new grid
    GenerateGrid();
    // generate new particles
    GenerateParticles();
    // place them in a gridcell
    AssignCellParticles();
    // and paint the new scene
    update();
}

// step through the simulation
void SimulationWidget::stepSimulation() {
    // step once
    LeapFrogIntegrate();
    // update cells
    AssignCellParticles();
    // paint result
    update();
}









