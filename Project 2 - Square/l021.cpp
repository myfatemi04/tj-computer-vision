#include<iostream>
#include<fstream>
#include<time.h>
#include<math.h>
#include<vector>
#include "geometry.cpp"
#include "image.cpp"
#include "color.cpp"

using namespace std;

void part2() {
    Polygon myQuad = Polygon::generateConvex(4);
    // Polygon myQuad = Polygon::fromFile("points.txt", 4);
    int width = 200;
    Image outputImage(width, width, 3, Colors::WHITE);

    /*
           B
    A--__ |
         X--C
     D   |
         E

    1) Drop a line BE perpendicular to AC such that BE = AC
    */
    Point A = myQuad.getPoint(0);
    Point B = myQuad.getPoint(1);
    Point C = myQuad.getPoint(2);
    Point D = myQuad.getPoint(3);
    LineSegment AC = LineSegment(A, C);
    LineSegment BE = AC.rotateClockwiseAroundA() + (B - A);
    Point E = BE.getB();

    outputImage.drawLineSegment(AC * width, Colors::BLUE);
    outputImage.drawLineSegment(BE * width, Colors::GREEN);

    Line sideDE = Line::fromPoints(D, E);
    Line sideA = sideDE.getPerpendicular().through(A);
    Line sideC = sideDE.getPerpendicular().through(C);
    Line sideB = sideDE.through(B);

    std::vector<Point> squarePoints;
    squarePoints.push_back(sideA.intersection(sideB));
    squarePoints.push_back(sideB.intersection(sideC));
    squarePoints.push_back(sideC.intersection(sideDE));
    squarePoints.push_back(sideDE.intersection(sideA));

    std::vector<Line> squareLines;

    for (int i = 0; i < 4; i++) {
        Point a = squarePoints.at(i);
        Point b = squarePoints.at((i + 1) % 4);
        squareLines.push_back(LineSegment(a, b).getLine());
        outputImage.drawCircle(Circle(a * width, 2), Colors::RED);
        outputImage.drawLine(LineSegment(a * width, b * width).getLine(), Colors::RED);
        // outputImage.drawLineSegment(a * width, b * width, Colors::RED);
        cout << a.getX() << ", " << a.getY() << "\n";
    }

    Polygon square(squarePoints);

    // Draw the quadrilateral
    outputImage.drawPolygon(myQuad * width, Colors::BLACK, true);
    
    ofstream outfile("output.ppm");

    outputImage.saveImage(outfile);

    outfile.close();
    
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

void part1() {
    Polygon myConvexQuad = Polygon::generateConvex(4);
    myConvexQuad.save("points.txt");
}

int main() {
    srand(time(nullptr));

    // part1();
    part2();
}