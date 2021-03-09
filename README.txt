setup of the application:
using qt and opengl to draw particles on to.... using the draw pixels method ?
an openglwidget
the simulation uses a 2D grid, so maybe use the textures approach
a gui specifying the initial state of the simulation (water and environment) ... viscosity, pressure, gravity

compute only with gravity and internal pressure at first.

Leapfrog scheme
We compute the acceleration with the forces.
v(n+1/2) = v(n) + a(n) * deltaT/2
x(n+1) = x(n) + v(n+1/2)
v(n+1) = v(n+1/2) + a(n+1) * deltaT/2

Divide the space up into a grid... we need a way to transform from grid to texture.
Maybe not directly influence h, but compute it by dividing a certain area (lxl) into a certain 
number of cells m such that h = l / m

Texture bottom left can be the origin (0,0)

At each iteration, we need to assign n particles to their corresponding cell O(n). 
Then we have iterate over all the cells and apply the search only to cell and its neighbours.