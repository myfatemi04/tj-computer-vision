#include<iostream>
#include<fstream>
#include<time.h>
#include<math.h>
#include<vector>
#include "geometry.cpp"
#include "image.cpp"
#include "color.cpp"

using namespace std;

void part2(Polygon myQuad) {
    int width = 800;
    Image outputImage(width, width, 3, Colors::WHITE);

    /*
           B
    A--__ |
         X--C
     D   |
         E

    1) Drop a line BE perpendicular to AC such that BE = AC
    */
    std::vector<Polygon> squares;
    std::vector<Line> smallestSquareLines;
    Polygon* smallestSquare = 0;
    Polygon* square;

    for (int iter = 0; iter < 6; iter++) {
        Point A = myQuad.getPoint(0);
        std::vector<int> available = {1, 2, 3};
        Point opp = myQuad.getPoint(available.at(iter % 3));

        // pop iter % 3rd
        available.erase(available.begin() + (iter % 3));

        // pop first
        Point other = myQuad.getPoint(available.at(0));
        available.erase(available.begin());

        // pop first
        Point last = myQuad.getPoint(available.at(0));
        available.erase(available.begin());

        LineSegment AC = LineSegment(A, opp);
        LineSegment BE = ((iter % 2 == 0) ? AC.rotateClockwiseAroundA() : AC.rotateCounterClockwiseAroundA()) + (other - A);
        Point E = BE.getB();

        Line sideDE = Line::fromPoints(last, E);
        Line sideA = sideDE.getPerpendicular().through(A);
        Line sideC = sideDE.getPerpendicular().through(opp);
        Line sideB = sideDE.through(other);

        std::vector<Point> squarePoints;
        squarePoints.push_back(sideA.intersection(sideB));
        squarePoints.push_back(sideB.intersection(sideC));
        squarePoints.push_back(sideC.intersection(sideDE));
        squarePoints.push_back(sideDE.intersection(sideA));

        // outputImage.drawLine(sideA * width, Colors::BLACK);
        // outputImage.drawLine(sideB * width, Colors::BLACK);
        // outputImage.drawLine(sideC * width, Colors::BLACK);
        // outputImage.drawLine(sideDE * width, Colors::BLACK);

        square = new Polygon(squarePoints);
        squares.push_back(*square);

        if (smallestSquare == NULL) {
            // std::cout << "Set smallest square on iter " << iter << "\n";
            smallestSquare = square;
            smallestSquareLines = {
                sideA, sideB, sideC, sideDE
            };
        } else {
            if ((*square).getSide(0).length() < (*smallestSquare).getSide(0).length()) {
                // std::cout << "Set smallest square on iter " << iter << "\n";
                smallestSquare = square;
                smallestSquareLines = {
                    sideA, sideB, sideC, sideDE
                };
            }
        }
    }

    for (int i = 0; i < 4; i++) {
        outputImage.drawLine(smallestSquareLines.at(i) * width, Colors::BLACK);
        outputImage.drawCircle(Circle(myQuad.getPoint(i) * width, 2), Colors::BLACK);
    }
    outputImage.drawPolygon((*smallestSquare) * width, Colors::PURPLE, 5);
    
    ofstream outputImageFile("output.ppm");
    outputImage.saveImage(outputImageFile);
    outputImageFile.close();

    ofstream outputPointsFile("output.txt");
    outputPointsFile.precision(17);
    outputPointsFile << myQuad << "\n";

    for (int i = 0; i < squares.size(); i++) {
        double squareSideLength = squares.at(i).getSide(0).length();
        double squareArea = squareSideLength * squareSideLength;
        outputPointsFile << squares.at(i) << " Area = " << squareArea << "\n";
    }

    outputPointsFile.close();
}

void part1() {
    Polygon myConvexQuad = Polygon::generateConvex(4);
    ofstream outputPointsFile("points.txt");
    outputPointsFile << myConvexQuad;
    outputPointsFile.close();
}

int main() {
    srand(time(nullptr));

    // part1();
    part2(Polygon::fromFile("points.txt", 4));
    // part2(Polygon::generateConvex(4));
}