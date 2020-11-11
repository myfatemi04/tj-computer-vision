#include <algorithm>
#include <iomanip>
#include <vector>
#include <cstring>
#include "geometry.cpp"
#include "image.cpp"

/// NUMBER OF POINTS TO USE
int npoints = 50;

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
 * General method to read `npoints` points from the input stream `in`.
 */
std::vector<Point> readPoints(std::istream& in, int npoints) {
  std::vector<Point> points;
  for (int i = 0; i < npoints; i++) {
    double x, y;
    in >> x >> y;
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
 * General method to save a vector of points to an output stream.
 */
void savePoints(std::ostream& out, std::vector<Point> points) {
  // makes all numbers write with exactly 17 digits after the decimal point
  out << std::fixed << std::setprecision(17);
  for (auto it = points.begin(); it != points.end(); it++) {
    out << it->getX() << "  " << it->getY() << "\n";
  }
}

void generatePoints() {
  std::ofstream pointsFile("points.txt");
  savePoints(pointsFile, generatePoints(npoints));
  pointsFile.close();
}

std::vector<Point> readPoints() {
  std::ifstream pointsFile("points.txt");
  std::vector<Point> points = readPoints(pointsFile, npoints);
  pointsFile.close();
  return points;
}

PointPair part1(std::vector<Point>& points) {
  PointPair closest(points.at(0), points.at(1));

  for (auto i = points.begin(); i != points.end(); i++) {
    for (auto j = std::next(i); j != points.end(); j++)  {
      closest.minify({*i, *j});
    }
  }

  return closest;
}

/**
 * Method to use when sorting points by X value
 */
bool comparePointXValues(const Point& first, const Point& second) {
  return first.getX() < second.getX();
}

/**
 * Divide and conquer method:
 *  finds the closest pair of points with both in the first half,
 *  finds the closest pair of points with both in the second half,
 *  and checks to see if the closest pair of points is shared between
 *  the first and second halves.
 */
PointPair helper(std::vector<Point>& points, int begin, int end) {
  int numPoints = (end - begin);
  int mid = (begin + end) >> 1;
  if (numPoints == 2) {
    return PointPair(points.at(begin), points.at(end - 1));
  } else if (numPoints == 3) {
    PointPair p1 = PointPair(points.at(begin + 0), points.at(begin + 1));
    PointPair p2 = PointPair(points.at(begin + 1), points.at(begin + 2));
    PointPair p3 = PointPair(points.at(begin + 2), points.at(begin + 0));
    
    PointPair res = p1;
    res.minify(p2);
    res.minify(p3);

    return res;
  } else {
    // Right side
    PointPair min1 = helper(points, begin, mid);
    // Left side
    PointPair min2 = helper(points, mid, end);
    // Closest pair where both are in left, or both are in right
    PointPair closest = min1;
    closest.minify(min2);

    // Maybe one of the points is on the left side and one is on the right?
    // Make a 'strip' of points down the middle
    int stripLeft = mid - 1, stripRight = mid;
    double middleX = points.at(mid).getX();

    // [begin, mid) = left side
    for (int i = mid - 1; i >= begin; i--) {
      double distanceFromMiddle = abs(middleX - points.at(i).getX());
      if (distanceFromMiddle <= closest.getDistance()) {
        stripLeft = i;
      } else {
        break; 
      }
    }

    // [mid, end] = right side
    for (int i = mid; i < end; i++) {
      double distanceFromMiddle = abs(middleX - points.at(i).getX());
      if (distanceFromMiddle <= closest.getDistance()) {
        stripRight = i;
      } else {
        break; 
      }
    }

    // Brute force the strip
    for (int leftPoint = stripLeft; leftPoint < mid; leftPoint++) {
      for (int rightPoint = mid; rightPoint < stripRight; rightPoint++) {
        closest.minify({
          points.at(leftPoint),
          points.at(rightPoint)
        });
      }
    }

    return closest;
  }
}

/**
 * Driver method for the Divide and conquer method:
 *  sorts the input vector, calls the helper method with the
 *  appropriate values
 */
PointPair part2(std::vector<Point>& points) {
  // First, sort the points by X value
  std::sort(points.begin(), points.end(), comparePointXValues);

  // Now, divide-and-conquer
  return helper(points, 0, points.size());
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

/**
 * Given `points`, find the closest pair of points using brute force.
 * Time how long it takes.
 */
void timePart1(std::vector<Point>& points, std::ofstream& outfile) {
  auto start = getMillis();
  PointPair closest = part1(points);
  auto elapsed = getMillis() - start;

  std::cout << "Brute force: " << elapsed << "ms\n";
  std::cout << closest << std::endl;
  
  outfile << "Brute force: " << elapsed << "ms\n";
  outfile << closest << std::endl;
}

/**
 * Given `points`, find the closest pair of points using divide and conquer.
 * Time how long it takes.
 */
void timePart2(std::vector<Point>& points, std::ofstream& outfile) {
  auto start = getMillis();
  PointPair closest = part2(points);
  auto elapsed = getMillis() - start;

  std::cout << "Recursive: " << elapsed << "ms\n";
  std::cout << closest << std::endl;

  outfile << "Recursive: " << elapsed << "ms\n";
  outfile << closest << std::endl;
}

/**
 * Main method. Sets the random seed to the current time.
 * Also takes in an argument for the number of points to use.
 * If not provided, the default number of points is 10000.
 * Reads the list of points from the default file.
 */
int main(int argc, const char* argv[]) {
  std::srand(time(NULL));
  if (argc > 1) {
    npoints = atoi(argv[1]);
  }
  
  // generatePoints();

  auto points = readPoints();

  std::ofstream outfile("results.txt");
  timePart1(points, outfile);
  timePart2(points, outfile);
  outfile.close();
}