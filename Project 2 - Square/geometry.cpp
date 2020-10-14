#include<math.h>
#include<iostream>
#include<fstream>

bool between(double a, double b, double p) {
    if (a < b) {
        return a <= p && p <= b;
    } else {
        return b <= p && p <= a;
    }
}

double getRandom() {
    return rand() / ((double) RAND_MAX);
}

double* pointFromAngle(double angle, double mag) {
    return new double[2] {mag * cos(angle), mag * sin(angle)};
}

class Point {
    private:
        double x, y;

    public:
        Point(double x, double y) {
            this -> x = x;
            this -> y = y;
        }

        Point() {}

        static double distance(Point a, Point b) {
            return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
        }

        double magnitude() {
            return sqrt(x * x + y * y);
        }

        Point operator+(Point b) {
            return Point(x + b.x, y + b.y);
        }

        Point operator-(Point b) {
            return Point(x - b.x, y - b.y);
        }

        Point operator*(double scalar) {
            return Point(x * scalar, y * scalar);
        }

        double getX() {
            return x;
        }

        double getY() {
            return y;
        }
};

void printPoint(FILE* fptr, Point point) {
    fprintf(fptr, "(%.17f,%.17f)", point.getX(), point.getY());
}

class LineSegment {
    private:
        Point a, b;

    public:
        LineSegment(Point a, Point b) {
            this -> a = a;
            this -> b = b;
        }

        LineSegment(double x1, double y1, double x2, double y2) {
            this -> a = Point(x1, y1);
            this -> b = Point(x2, y2);
        }

        double length() {
            return Point::distance(a, b);
        }

        double slope() {
            return (b.getY() - a.getY()) / (b.getX() - a.getX());
        }

        bool pointAbove(Point p) {
            // vertical line must intersect
            // first, finds the slope of the line
            if (!between(a.getX(), b.getX(), p.getX())) {
                return false;
            }

            // equation: y = mx + b
            // a[1] = ma[0] + b
            double intercept = a.getY() - slope() * a.getX();
            double segmentY = slope() * p.getX() + intercept;
            return p.getY() > segmentY;
        }

        LineSegment extendTo(double length) {
            Point offset = b - a;
            offset = offset * (length / this -> length());
            return LineSegment(a, a + offset);
        }
};

class Triangle {
    private:
        Point points[3];
        LineSegment* sides;

    public:
        Triangle(Point* points) {
            this -> points[0] = points[0];
            this -> points[1] = points[1];
            this -> points[2] = points[2];

            sides = new LineSegment[3] {
                LineSegment(points[0], points[1]),
                LineSegment(points[1], points[2]),
                LineSegment(points[2], points[0])
            };
        }

        LineSegment* getSides() {
            return sides;
        }

        bool containsPoint(Point p) {
            int count = 0;
            for (int i = 0; i < 3; i++) {
                if (sides[i].pointAbove(p)) {
                    count += 1;
                }
            }

            return (count % 2) != 0;
        }
};

class Quadrilateral {
    public:
        Point* points;

        Quadrilateral(Point* points) {
            for (int i = 0; i < 4; i++) {
                this -> points[i] = points[i];
            }
        }

        bool isConvex() {
            for (int checkPointIndex = 0; checkPointIndex < 4; checkPointIndex++) {
                Point* trianglePoints = new Point[3];
                Point checkPoint = points[checkPointIndex];
                int pointsRecorded = 0;

                for (int trianglePointIndex = 0; trianglePointIndex < 4; trianglePointIndex++) {
                    if (trianglePointIndex != checkPointIndex) {
                        trianglePoints[pointsRecorded] = points[trianglePointIndex];
                        pointsRecorded += 1;
                    }
                }

                Triangle t = Triangle(trianglePoints);

                if (t.containsPoint(checkPoint)) {
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
            std::ifstream file(filename);
            Point* points = new Point[4];
            for (int i = 0; i < 4; i++) {
                double x, y;
                file.ignore(1);
                file >> x;
                file.ignore(1);
                file >> y;
                file.ignore(4);
                points[i] = Point(x, y);
            }

            file.close();

            return Quadrilateral { points };
        }

        /**
         * Guaranteed to generate in clockwise order
         * 1) Generates angles in clockwise order
         * 2) Generates magnitudes randomly
         * 3) Generates points based on angles and magnitudes
         */
        static Quadrilateral getRandomQuadrilateral() {
            double* slices = new double[4];
            double sum = 0;
            double circleRads = 6.28318531;

            for (int i = 0; i < 4; i++) {
                // generate proportions for 'slices'
                slices[i] = getRandom();

                // keep track of the sum of the slices so we can scale them to 2pi
                sum += slices[i];
            }
            
            Point* points = new Point[4];
            double angle = getRandom();
            for (int i = 0; i < 4; i++) {
                angle += circleRads * slices[i] / sum;

                // scale cos/sin to reach the square
                double maxMagnitude;
                double width = abs(cos(angle));
                double height = abs(sin(angle));
                if (width > height) {
                    maxMagnitude = 1 / width;
                } else {
                    maxMagnitude = 1 / height;
                }

                double magnitude = maxMagnitude * getRandom();

                // Generate with magnitude/direction
                points[i] = Point(magnitude * cos(angle), magnitude * sin(angle));
                
                // Fit to the unit square
                points[i] = (points[i] * 0.5) + Point(0.5, 0.5);

                std::cout << "Generating point; m=" << magnitude;
                std::cout << "; theta=" << (angle * 180 / 3.14159265) << "\n";
            }

            return points;
        }

};
