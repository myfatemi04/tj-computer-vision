#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>

#ifndef GEOMETRY
#define GEOMETRY

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
        Point(double x, double y): x(x), y(y) {}

        Point(const Point& point): x(point.x), y(point.y) {}

        Point() {}

        static double distance(const Point& a, const Point& b) {
            return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
        }

        double distance(const Point& other) const {
            return Point::distance(*this, other);
        }

        double magnitude() const {
            return sqrt(x * x + y * y);
        }

        Point operator+(const Point& b) const {
            return Point(x + b.x, y + b.y);
        }

        Point operator-(const Point& b) const {
            return Point(x - b.x, y - b.y);
        }

        Point operator*(double scalar) const {
            return Point(x * scalar, y * scalar);
        }

        double getX() const { return x; }
        double getY() const { return y; }
};

void printPoint(FILE* fptr, Point point) {
    fprintf(fptr, "(%.17f,%.17f)", point.getX(), point.getY());
}

class Line {
    private:
        double a, b, c; // ax + by = c

    public:
        Line() {}
        Line(
            double a,
            double b,
            double c
        ):
            a(a),
            b(b),
            c(c) {}

        double getA() const { return a; }
        double getB() const { return b; }
        double getC() const { return c; }

        double getSlope() const {
            return -a / b;
        }

        double getIntercept() const {
            return c / b;
        }

        bool isParallel(const Line& other) const {
            if (b == 0) {
                return other.b == 0;
            } else {
                if (other.b == 0) {
                    return false;
                } else {
                    return other.getSlope() == getSlope();
                }
            }
        }

        bool isPerpendicular(const Line& other) const {
            if (b == 0) {
                return other.a == 0;
            } else {
                if (other.a == 0) {
                    return false;
                } else {
                    return (other.b / other.a) == (-a / b);
                }
            }
        }

        Point intersection(const Line& other) const {
            // matrix is:
            /*
            | a1 b1 | = | c1 |
            | a2 b2 |   | c2 |

            We will use Cramer's rule to find the (x, y) of the intersection.
            */

            double det = a * other.b - b * other.a;

            if (det == 0) {
                throw 0;
            }

            double detX = c * other.b - b * other.c;
            double detY = a * other.c - c * other.a;

            return Point(detX / det, detY / det);
        }

        Line through(const Point& p) const {
            double targetC = (p.getX() * a) + (p.getY() * b);
            return Line(a, b, targetC);
        }

        Line getPerpendicular() const {
            // ax + by = c
            // by = -ax + c
            // ay = bx + c
            // -bx + ay = c
            return Line(-b, a, c);
        }

        static Line fromPoints(Point a, Point b) {
            
            /*
            equations are ax + by = c
            how to find a, b, and c?
            y = mx + b
            m = slope_numer/slope_denom
            y = slope_numer/slope_denom * x + b
            slope_denom(y) = slope_numer(x) + b
            -slope_numer(x) + slope_denom(y) = b
            now, substitute a point for (x, y) to find b. we can just use (x1, y1)
            b = -slope_numer(x1) + slope_denom(y1)
            now, we know that in ax + by = c:
            a = -slope_numer
            b = slope_denom
            c = a(x1) + b(y1)
            */
            double slopeNumer, slopeDenom;
            slopeNumer = b.getY() - a.getY();
            slopeDenom = b.getX() - a.getX();

            double newA = -slopeNumer;
            double newB = slopeDenom;
            double newC = newA * (a.getX()) + newB * (a.getY());

            return Line(newA, newB, newC);
        }

        Line operator*(double scalar) const {
            return Line(a, b, c * scalar);
        }
};

class LineSegment {
    private:
        Point a, b;

    public:
        LineSegment(Point a, Point b): a(a), b(b) {}

        LineSegment(double x1, double y1, double x2, double y2):
            a(Point(x1, y2)),
            b(Point(x2, y2)) {}

        double length() const {
            return Point::distance(a, b);
        }

        double slope() const {
            return (b.getY() - a.getY()) / (b.getX() - a.getX());
        }

        bool pointAbove(Point p) const {
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

        Line getLine() const {
            return Line::fromPoints(a, b);
        }

        LineSegment rotateClockwiseAroundA() {
            Point offset = b - a;
            Point newOffset = Point(offset.getY(), -offset.getX());
            return LineSegment(a, a + newOffset);
        }

        LineSegment rotateClockwiseAroundB() {
            Point offset = a - b;
            Point newOffset = Point(offset.getY(), -offset.getX());
            return LineSegment(b + newOffset, b);
        }

        LineSegment rotateCounterClockwiseAroundA() {
            Point offset = b - a;
            Point newOffset = Point(-offset.getY(), offset.getX());
            return LineSegment(a, a + newOffset);
        }

        LineSegment rotateCounterClockwiseAroundB() {
            Point offset = a - b;
            Point newOffset = Point(-offset.getY(), offset.getX());
            return LineSegment(b + newOffset, b);
        }

        LineSegment operator+(Point p) {
            return LineSegment(a + p, b + p);
        }

        LineSegment operator-(Point p) {
            return LineSegment(a - p, b - p);
        }

        LineSegment operator*(double scalar) {
            return LineSegment(a * scalar, b * scalar);
        }

        Point getA() const { return a; }
        Point getB() const { return b; }
};

class Polygon {
    public:
        std::vector<Point> points;
        std::vector<LineSegment> sides;

        Polygon(std::vector<Point> points) {
            for (int i = 0; i < points.size(); i++) {
                this -> points.push_back(points.at(i));
                this -> sides.push_back(LineSegment(points.at(i), points.at((i + 1) % points.size())));
            }
        }

        Point getPoint(int point) const {
            return points.at(point);
        }

        LineSegment getSide(int side) const {
            return sides.at(side);
        }

        std::vector<Point> getPoints() const {
            return points;
        }

        std::vector<LineSegment> getSides() const {
            return sides;
        }

        bool containsPoint(Point p) const {
            int count = 0;
            for (int i = 0; i < (sides.size() - 1); i++) {
                if (sides.at(i).pointAbove(p)) {
                    count += 1;
                }
            }

            return (count % 2) != 0;
        }

        bool isConvex() const {
            for (int checkPointIndex = 0; checkPointIndex < 4; checkPointIndex++) {
                std::vector<Point> otherPoints;
                Point checkPoint = points.at(checkPointIndex);
                for (int otherPointIndex = 0; otherPointIndex < 4; otherPointIndex++) {
                    if (otherPointIndex != checkPointIndex) {
                        otherPoints.push_back(points[otherPointIndex]);
                    }
                }

                Polygon otherPointsPolygon(otherPoints);

                if (otherPointsPolygon.containsPoint(checkPoint)) {
                    return false;
                }
            }

            return true;
        }

        static Polygon generateConvex(int nsides) {
            while (true) {
                std::cout << "Generating...\n";
                Polygon tmp (getRandomPolygon(nsides));

                if (tmp.isConvex()) {
                    return tmp;
                }
            }
        }

        static Polygon getRandomPolygon(int nsides) {
            std::vector<Point> points;
            for (int i = 0; i < nsides; i++) {
                points.push_back(Point(getRandom(), getRandom()));
            }
            return Polygon(points);
        }

        /**
         * Guaranteed to generate in counterclockwise order
         * 1) Generates angles in counterclockwise order
         * 2) Generates magnitudes randomly
         * 3) Generates points based on angles and magnitudes
         */
        static Polygon getRandomPolygonCounterclockwise(int nsides) {
            double* slices = new double[nsides];
            double sum = 0;
            double circleRads = 6.28318531;

            for (int i = 0; i < nsides; i++) {
                // generate proportions for 'slices'
                slices[i] = getRandom();

                // keep track of the sum of the slices so we can scale them to 2pi
                sum += slices[i];
            }
            
            std::vector<Point> points;
            double angle = getRandom();
            for (int i = 0; i < nsides; i++) {
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
                Point newPoint(magnitude * cos(angle), magnitude * sin(angle));
                // Fit to the unit square
                newPoint = (newPoint * 0.5) + Point(0.5, 0.5);

                points.push_back(newPoint);

                std::cout << "Generating point; m=" << magnitude;
                std::cout << "; theta=" << (angle * 180 / 3.14159265) << "\n";
            }

            return Polygon(points);
        }

        static Polygon fromFile(const char* filename, int nsides) {
            std::ifstream file(filename);
            std::vector<Point> points;
            for (int i = 0; i < nsides; i++) {
                double x, y;
                file.ignore(1); // '('
                file >> x;
                file.ignore(1); // ','
                file >> y;
                file.ignore(1); // ')'

                if (i < nsides - 1) {
                    file.ignore(3); // ' , '
                }

                points.push_back(Point(x, y));
            }

            file.close();

            return Polygon(points);
        }

        int length() const {
            return points.size();
        }

        double perimeter() const {
            double acc = 0;
            for (auto it : sides) {
                acc += it.length();
            }
            return acc;
        }

        Polygon operator*(double scalar) {
            std::vector<Point> newPoints;
            for (int i = 0; i < points.size(); i++) {
                newPoints.push_back(getPoint(i) * scalar);
            }
            return Polygon(newPoints);
        }
};

class Circle {
    private:
        Point center;
        double radius;

    public:
        Circle(Point center, double radius) {
            this -> center = center;
            this -> radius = radius;
        }

        double getRadius() const {
            return this -> radius;
        }

        Point getCenter() const {
            return this -> center;
        }

        Circle operator*(double scalar) {
            return Circle(center * scalar, radius * scalar);
        }
};

// CUSTOM OUTPUT METHODS
// POINT
std::ostream& operator<<(std::ostream &stream, const Point &point) {
    return stream << "(" << point.getX() << "," << point.getY() << ")";
}

// LINE
std::ostream& operator<<(std::ostream &stream, const Line &line) {
    return stream << "<" << line.getA() << ", " << line.getB() << ", " << line.getC() << ">";
}

// LINE SEGMENT
std::ostream& operator<<(std::ostream &stream, const LineSegment &lineSegment) {
    return stream << "<" << lineSegment.getA() << ", " << lineSegment.getB() << ">";
}

// POLYGON
std::ostream& operator<<(std::ostream &stream, const Polygon &polygon) {
    std::vector<Point> points = polygon.getPoints();
    for (int i = 0; i < polygon.length(); i++) {
        stream << points.at(i);
        if (i < polygon.length() - 1) {
            stream << " , ";
        }
    }

    return stream;
}

// CIRCLE

std::ostream& operator<<(std::ostream &stream, const Circle &circle) {
    stream << "Circle" << circle.getCenter() << " r=" << circle.getRadius();
    return stream;
}

#endif