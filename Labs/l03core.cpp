#ifndef LAB_03_CORE
#define LAB_03_CORE

#include <algorithm>
#include <iomanip>
#include <vector>
#include "geometry.cpp"

/**
 * A class that stores a pair of points and their distance.
 * Includes a method to find the minimum and replace the points
 * with the argument if the argument's pair of points are closer
 * together. Utilizes references to spend less time copying by value.
 * Michael Fatemi 11/6/2020
 */
class PointPair {
  private:
    Point *a, *b;
    double distance;

  public:
    PointPair(Point& a, Point& b): a(&a), b(&b) {
      distance = Point::distance(a, b);
    }

    double getDistance() const {
      return distance;
    }

    Point getA() const { return *a; }
    Point getB() const { return *b; }
    
    void minify(PointPair other) {
      if (other.getDistance() < distance) {
        a = other.a;
        b = other.b; 
        distance = other.distance;
      }
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
void savePoints(std::vector<Point>& points, const char* filename = "points.txt") {
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

#include <chrono>
/**
 * Uses the `chrono` module to get the system time in milliseconds.
 * Allows us to get the elapsed runtime of each method.
 */
long long getMillis() {
  std::chrono::milliseconds ms = std::chrono::duration_cast< std::chrono::milliseconds >(
    std::chrono::system_clock::now().time_since_epoch()
  );

  return ms.count();
}

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

  std::cout << methodName << " x [" << count << "] ";
  std::cout << elapsed << "ms (" << percycle << "/cycle)\n";
  std::cout << closest << std::endl;

  outfile << methodName << " x [" << count << "] ";
  outfile << elapsed << "ms (" << percycle << "/cycle)\n";
  outfile << closest << std::endl;
}

#endif