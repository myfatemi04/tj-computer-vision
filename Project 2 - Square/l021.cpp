#include<iostream>
#include<fstream>
#include<time.h>

using namespace std;

void printPoint(FILE* fptr, double* point) {
    fprintf(fptr, "(%.17f,%.17f)", point[0], point[1]);
}

bool pointInsideTriangle(double** triangle, double* pointToCheck) {
    return false;
}

double getRandom() {
    return rand() / ((double) RAND_MAX);
}

double** getRandomQuadrilateral() {
    double** points = new double*[4];
    for (int i = 0; i < 4; i++) {
        points[i] = new double[2];
        points[i][0] = getRandom();
        points[i][1] = getRandom();
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

int main() {
    srand(time(nullptr));

    Quadrilateral q = Quadrilateral::fromFile("points.txt");
    q.save("points_copy.txt");
    // Quadrilateral myConvexQuad = Quadrilateral::generateConvex();
    // myConvexQuad.save("points.txt");
}