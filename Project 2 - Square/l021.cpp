#include<iostream>
#include<fstream>
#include<time.h>
#include "geometry.cpp"

using namespace std;

void part1() {
    Polygon myConvexQuad = Polygon::generateConvex(4);
    ofstream outputPointsFile("points.txt");
    outputPointsFile << myConvexQuad;
    outputPointsFile.close();
}

int main() {
    srand(time(nullptr));

    part1();
}