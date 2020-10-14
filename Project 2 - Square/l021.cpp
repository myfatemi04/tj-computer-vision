#include<iostream>
#include<fstream>
#include<time.h>
#include<math.h>
#include "image.cpp"

using namespace std;

void printPoint(FILE* fptr, double* point) {
    fprintf(fptr, "(%.17f,%.17f)", point[0], point[1]);
}

bool between(double a, double b, double p) {
    if (a < b) {
        return a <= p && p <= b;
    } else {
        return b <= p && p <= a;
    }
}

bool pointAbove(double* a, double* b, double* p) {
    // vertical line must intersect
    // first, finds the slope of the line
    if (between(a[0], b[0], p[0])) {
        // equation: y = mx + b
        // a[1] = ma[0] + b
        double slope = (b[1] - a[1]) / (b[0] - a[0]);
        double intercept = a[1] - slope * a[0];
        double segmentY = slope * p[0] + intercept;
        return p[1] > segmentY;
    } else {
        return false;
    }
}

bool pointInsideTriangle(double** triangle, double* pointToCheck) {
    int count = 0;
    if (pointAbove(triangle[0], triangle[1], pointToCheck)) count += 1;
    if (pointAbove(triangle[1], triangle[2], pointToCheck)) count += 1;
    if (pointAbove(triangle[2], triangle[0], pointToCheck)) count += 1;

    if (count == 0 || count == 2) {
        return false;
    } else {
        return true;
    }
}

double getRandom() {
    return rand() / ((double) RAND_MAX);
}

double* pointFromAngle(double angle, double mag) {
    return new double[2] {mag * cos(angle), mag * sin(angle)};
}

/**
 * Guaranteed to generate in clockwise order
 * 1) Generates angles in clockwise order
 * 2) Generates magnitudes randomly
 * 3) Generates points based on angles and magnitudes
 */
double** getRandomQuadrilateral() {
    double* slices = new double[4];
    double sum = 0;
    double circleRads = 6.28318531;

    for (int i = 0; i < 4; i++) {
        // generate proportions for 'slices'
        slices[i] = getRandom();

        // keep track of the sum of the slices so we can scale them to 2pi
        sum += slices[i];
    }

    double** points = new double*[4];
    double angle = getRandom();
    for (int i = 0; i < 4; i++) {
        angle += circleRads * slices[i] / sum;

        // scale cos/sin to reach the square
        double maxMagnitude;
        if (cos(angle) > sin(angle)) {
            maxMagnitude = 1 / cos(angle);
        } else {
            maxMagnitude = 1 / sin(angle);
        }

        double magnitude = maxMagnitude * getRandom();

        // Generate with magnitude/direction
        points[i] = new double[2] {magnitude * cos(angle), magnitude * sin(angle)};

        // Fit to the unit square
        points[i][0] = points[i][0] * 0.5 + 0.5;
        points[i][1] = points[i][1] * 0.5 + 0.5;
    }

    return points;
}

class Quadrilateral {
    public:
        double** points;

        Quadrilateral(double** points_): points(points_) {}

        bool isConvex() {
            for (int checkPointIndex = 0; checkPointIndex < 4; checkPointIndex++) {
                double** trianglePoints = new double*[3];
                double* checkPoint = points[checkPointIndex];
                int pointsRecorded = 0;

                for (int trianglePointIndex = 0; trianglePointIndex < 4; trianglePointIndex++) {
                    if (trianglePointIndex != checkPointIndex) {
                        trianglePoints[pointsRecorded] = points[trianglePointIndex];
                        pointsRecorded += 1;
                    }
                }

                if (pointInsideTriangle(trianglePoints, checkPoint)) {
                    return false;
                }
            }

            return true;
        }

        static Quadrilateral generateConvex() {
            while (true) {
                Quadrilateral tmp (getRandomQuadrilateral());

                if (tmp.isConvex()) {
                    return tmp;
                }
            }
        }

        void save(const char* filename) {
            FILE* fptr = fopen(filename, "w");
            for (int i = 0; i < 4; i++) {
                printPoint(fptr, points[i]);
                if (i < 3) {
                    fprintf(fptr, " , ");
                }
            }
        }

        static Quadrilateral fromFile(const char* filename) {
            ifstream file(filename);
            double** points = new double*[4];
            for (int i = 0; i < 4; i++) {
                points[i] = new double[2];
                file.ignore(1);
                file >> points[i][0];
                file.ignore(1);
                file >> points[i][1];
                file.ignore(4);

                cout << points[i][0] << ", " << points[i][1] << endl;
            }

            file.close();

            return Quadrilateral { points };
        }
};

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