#include<iostream>
#include<fstream>
#include<time.h>
#include<math.h>
#include "geometry.cpp"
#include "image.cpp"

using namespace std;

void part2() {
    // Quadrilateral myQuad = Quadrilateral::fromFile("points.txt");
    // int white[3] = {1, 1, 1};
    // Image outputImage(800, 800, 3, white);

    // Create all possible squares

    
    // Display all 4 points of the quad by creating a circle with radius=2
    // Draw square with the minimum area
    // Save the coordinates of the 4 points as well as 4 vertices of all squares found
    // Use the same precision as with the square
    /*

    (p1,p1) , (p2,p2) , (p3,p3) , (p4,p4)
    (vx1,vy1) , (vx2,vy2) , (vx3,vy3) , (vx4,vy4) Area=ar1
    ....

    */
}

int main() {
    srand(time(nullptr));

    // Quadrilateral q = Quadrilateral::fromFile("points.txt");
    // q.save("points_copy.txt");
    Quadrilateral myConvexQuad = Quadrilateral::generateConvex();
    myConvexQuad.save("points.txt");
}