COMP 5823 Animation and Simulation cw4

Thomas Moreno Cooper sc17tm@leeds.ac.uk
University of Leeds 
mar 2021

Libraries used:
Qt

! uses c++ 11 !

With the libraries and Qt version 5.9.5 on Linux x86
- run qmake (version 3.1)
- make
- execute with ./Dungeon4

The simulation slows down for large numbers of particles as the application is not multi-threaded. Reducing the kernel radius for computing forces improves the computation time (iterate over fewer particles) at the cost of simulation stability. 

![Fluid sim](https://media.giphy.com/media/DiIu3nEeinEnp8YhYc/giphy.gif)
