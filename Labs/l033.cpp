#include <algorithm>
#include <vector>
#include "geometry.cpp"
#include "l03core.cpp"

/**
 * Improved divide and conquer method. Optimizes based on the fact
 *  that you can make the comparisons between strips work in O(n log n)
 *  time, which is an improvement over O(n) time.
 */
PointPair helper3(std::vector<Point>& points, int begin, int end) {
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
    PointPair min1 = helper3(points, begin, mid);
    // Left side
    PointPair min2 = helper3(points, mid, end);
    // Closest pair where both are in left, or both are in right
    PointPair closest = min1;
    closest.minify(min2);
    double d = closest.getDistance();

    // Maybe one of the points is on the left side and one is on the right?
    // Make a 'strip' of points down the middle
    int stripLeft = mid, stripRight = mid;
    double middleX = points.at(mid).getX();
    double minX = middleX - d;
    double maxX = middleX + d;

    // [begin, mid) = left side
    while ((stripLeft > begin) && (points.at(stripLeft - 1).getX() >= minX)) {
      stripLeft--;
    }

    // [mid, end] = right side
    while ((stripRight < end) && (points.at(stripRight).getX() <= maxX)) {
      stripRight++;
    }

    int stripSize = stripRight - stripLeft;

    // Create a strip sorted by y values
    std::vector<Point> strip(points.begin() + stripLeft, points.begin() + stripRight);
    std::sort(strip.begin(), strip.end(), comparePointYValues);

    // For each point in the strip, compare it with the next 15 points
    for (int thisPoint = 0; thisPoint < stripSize; thisPoint++) {
      for (
        int thatPoint = thisPoint + 1;
        (thatPoint < thisPoint + 8) && (thatPoint < stripSize);
        thatPoint++
      ) {
        closest.minify({strip.at(thisPoint), strip.at(thatPoint)});
      }
    }

    return closest;
  }
}

/**
 * Driver method for the Divide and conquer method:
 *  sorts the input vector, calls the helper method with the
 *  appropriate values. Optimized to only check the next 15 points.
 */
PointPair part3(std::vector<Point>& points) {
  // First, sort the points by X value
  std::sort(points.begin(), points.end(), comparePointXValues);

  // Now, divide-and-conquer
  return helper3(points, 0, points.size());
}

#ifndef LAB_03_MAIN
#define LAB_03_MAIN

/**
 * Main method. Sets the random seed to the current time.
 * Also takes in an argument for the number of points to use.
 * If not provided, the default number of points is 10000.
 * Reads the list of points from the default file.
 */
int main(int argc, const char* argv[]) {
  std::srand(time(NULL));

  if (argc > 1) {
    std::vector<Point> points = generatePoints(atoi(argv[1]));
    savePoints(points);
  }

  auto points = readPoints();

  std::ofstream outfile("results.txt");
  timer(points, outfile, part3, "Recursive Optimized", 10);
  outfile.close();
}

#endif