# Trailblazer
This project implements various algorithms to find paths on maps as well as generate random new ones. Please do not hesitate to contact me for any 
information regarding this on my email devansh@stanford.edu.
## Disclaimer
Stanford University's standard C++ libraries have been used to develop the project. I wrote the main implementation files i.e. 
[trailblazer.cpp] (https://github.com/devansh287/Trailblazer/blob/master/trailblazer.cpp) Trailblazer/src folder.
## Details (Taken from assignment instructions)
This program displays various 2-dimensional worlds that represent either maps, mazes, or terrain and allows the user to generate paths 
in a world from one point to another. 
If you click on any two points in the world, the program will find a path from the starting position to the ending position. 
As it does so, it will color the vertexes green, yellow, and gray based on the colors assigned to them by the algorithm. Once the 
path is found, the program will highlight it and display information about the path weight in the console. The user can select one of 
five path-searching algorithms in the top menu:

1. depth-first search (DFS)
2. breadth-first search (BFS)
3. Dijkstra's algorithm
4. A* search
5. Alternate Path
The window also contains several controls. You can load mazes and terrains of different sizes (tiny, small, medium, large, and huge) 
from the bottom drop-down menu and then clicking the "Load" button.

You can also generate random mazes using the Kruskal's algorithm.

