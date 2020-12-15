#ifndef LAB_PART
#define LAB_PART 2
#endif

#ifndef LAB_L032
#define LAB_L032

#include <algorithm>
#include <vector>
#include "l03core.cpp"

/**
 * Divide and conquer method:
 *  finds the closest pair of points with both in the first half,
 *  finds the closest pair of points with both in the second half,
 *  and checks to see if the closest pair of points is shared between
 *  the first and second halves.
 */
PointPair helper2(std::vector<Point>& points, int begin, int end) {
  int numPoints = (end - begin);
  int mid = (begin + end) >> 1;
  if (numPoints == 2) {
    return PointPair(points[begin], points[end - 1]);
  } else if (numPoints == 3) {
    PointPair closest;
    closest.minify(PointPair(points[begin + 0], points[begin + 1]));
    closest.minify(PointPair(points[begin + 1], points[begin + 2]));
    closest.minify(PointPair(points[begin + 2], points[begin + 0]));

    return closest;
  } else {
    PointPair closest;
    // Right side
    closest.minify(helper2(points, begin, mid));
    // Left side
    closest.minify(helper2(points, mid, end));
    // Closest pair where both are in left, or both are in right
    double d = closest.getDistance();

    // Maybe one of the points is on the left side and one is on the right?
    // Make a 'strip' of points down the middle
    int stripLeft = mid, stripRight = mid;
    double middleX = points[mid].getX();

    // [begin, mid) = left side
    while ((stripLeft > begin) && (points[stripLeft - 1].getX() >= middleX - d)) {
      stripLeft--;
    }

    // [mid, end] = right side
    while ((stripRight < end) && (points[stripRight].getX() <= middleX + d)) {
      stripRight++;
    }

    // For each point in the left side of the strip, compare it with
    // a point on the right side of the strip until the difference in
    // X is greater than d.
    for (int i = stripLeft; i < mid; i++) {
      double maxX = points[i].getX() + d;
      for (int j = mid; j < stripRight && points[j].getX() < maxX; j++) {
        closest.minify({
          points[i],
          points[j]
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
  return helper2(points, 0, points.size());
}

#endif

#if LAB_PART == 2

#include "l031.cpp"

/**
 * Main method. Sets the random seed to the current time.
 * Also takes in an argument for the number of points to use.
 * If not provided, the default number of points is 10000.
 * Reads the list of points from the default file.
 */
int main(int argc, const char* argv[]) {
  std::srand(time(NULL));

  std::vector<Point> points;

  if (argc > 1) {
    points = readPoints(argv[1]);
  } else {
    points = readPoints();
  }

  std::ofstream outfile("results.txt");
  // timer(points, outfile, part1, "Brute Force", 1);
  timer(points, outfile, part2, "Recursive", 10);
  outfile.close();
}

#endif