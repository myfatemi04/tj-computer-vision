#ifndef LAB_03_CORE
#define LAB_03_CORE

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>

double getRandom() {
    return rand() / ((double) RAND_MAX);
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

		/**
		 * Checks if two points are within a certain distance of each other, inclusive.
		 */
		static bool areWithinDistance(const Point& a, const Point& b, double distance) {
			return (
				((a.x - b.x) * (a.x - b.x)) + 
				((a.y - b.y) * (a.y - b.y))
			) <= (distance * distance);
		}

		double distance(const Point& other) const {
			return Point::distance(*this, other);
		}

		double magnitude() const {
			return sqrt(x * x + y * y);
		}

		double getX() const { return x; }
		double getY() const { return y; }
};

std::ostream& operator<<(std::ostream &stream, const Point &point) {
    return stream << "(" << point.getX() << "," << point.getY() << ")";
}

/**
 * A class that stores a pair of points and their distance.
 * Includes a method to find the minimum and replace the points
 * with the argument if the argument's pair of points are closer
 * together. Utilizes references to spend less time copying by value.
 * Michael Fatemi 11/6/2020
 */
class PointPair {
	private:
		Point a, b;
		double distance;
		bool initialized = false;

	public:
		PointPair() {}

		PointPair(Point a, Point b): a(a), b(b) {
			distance = Point::distance(a, b);
			initialized = true;
		}

		double getDistance() const {
			return distance;
		}

		Point getA() const { return a; }
		Point getB() const { return b; }
		
		bool minify(const PointPair& other) {
			if (
				other.initialized && (
					!initialized ||
					other.getDistance() < distance
				)
			) {
				a = other.a;
				b = other.b; 
				distance = other.distance;
				initialized = true;
				return true;
			}

			return false;
		}

		void minify(PointPair *other) {
			if (other != nullptr) {
				this -> minify(*other);
			}
		}

		bool isInitialized() const {
			return this -> initialized;
		}
};

/**
 * Make a PointPair writeable to an output stream.
 * Format: (point1) (point2) distance
 */
std::ostream& operator<<(std::ostream& s, const PointPair& p) {
	return s << p.getA() << " " << p.getB() << " " << p.getDistance();
}

/**
 * General method to save a vector of points to points.txt
 */
void savePoints(const std::vector<Point>& points, const char* filename = "points.txt") {
	std::ofstream out(filename);
	
	// makes all numbers write with exactly 17 digits after the decimal point
	out << std::fixed << std::setprecision(17);
	for (auto it = points.begin(); it != points.end(); it++) {
		out << it->getX() << "  " << it->getY() << "\n";
	}

	out.close();
}

/**
 * General method to read from points.txt until the end of the file
 */
std::vector<Point> readPoints(const char* filename = "points.txt") {
	std::ifstream in(filename);
	std::vector<Point> points;
	double x, y;
	while (in >> x) {
		in >> y;
		points.push_back(Point(x, y));
	}
	return points;
}

/**
 * General method to generate `npoints` points.
 */
std::vector<Point> generatePoints(int npoints) {
	std::vector<Point> points;
	for (int i = 0; i < npoints; i++) {
		points.push_back(Point(getRandom(), getRandom()));
	}
	return points;
}

/**
 * Method to use when sorting points by X value
 */
bool comparePointXValues(const Point& first, const Point& second) {
	return first.getX() < second.getX();
}

/**
 * Method to use when sorting points by Y value
 */
bool comparePointYValues(const Point& first, const Point& second) {
	return first.getY() < second.getY();
}

/**
 * Method to use when sorting points by Y value
 */
bool comparePointPointerYValues(const Point* first, const Point* second) {
	return first->getY() < second->getY();
}

#include <chrono>
/**
 * Uses the `chrono` module to get the system time in milliseconds.
 * Allows us to get the elapsed runtime of each method.
 */
long long getMillis() {
	std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::system_clock::now().time_since_epoch()
	);

	return ms.count();
}

const char* files[] = {
	"points_files/points100.txt",
	"points_files/points500.txt",
	"points_files/points1k.txt",
	"points_files/points2k.txt",
	"points_files/points5k.txt",
	"points_files/points10k.txt",
	"points_files/points20k.txt",
	"points_files/points50k.txt",
	"points_files/points100k.txt",
	"points_files/points200k.txt",
	"points_files/points500k.txt",
	"points_files/points1m.txt"
};

#include <functional>

/**
 * Timer method
 */
void timer(
	std::vector<Point>& points,
	std::ofstream& outfile,
	std::function<PointPair (std::vector<Point>&)> method,
	const char* methodName,
	int count = 1
) {
	auto start = getMillis();
	for (int i = 0; i < count - 1; i++) {
		method(points);
	}
	PointPair closest = method(points);
	auto elapsed = getMillis() - start;
	auto percycle = elapsed / count;

	// write to a stringstream so we can write the same
	// output to two other streams
	std::stringstream out;

	out << std::setprecision(17);
	out << methodName << " x [" << count << "] ";
	out << elapsed << "ms (" << percycle << "/cycle) " << points.size() << " points\n";
	out << closest << std::endl;

	std::cout << out.str();
	outfile << out.str();
}

#endif