#ifndef LAB_PART
#define LAB_PART 2
#endif

#ifndef LAB_L032
#define LAB_L032

#include <algorithm>
#include <vector>
#include "l03core.cpp"
#include "geometry.cpp"

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
    PointPair min1 = helper2(points, begin, mid);
    // Left side
    PointPair min2 = helper2(points, mid, end);
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

    // For each point in the left side of the strip, compare it with
    // a point on the right side of the strip until the difference in
    // X is greater than d.
    for (int i = stripLeft; i < mid; i++) {
      double maxX = points.at(i).getX() + d;
      for (int j = mid; j < stripRight && points.at(j).getX() < maxX; j++) {
        closest.minify({
          points.at(i),
          points.at(j)
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

/**
 * Main method. Sets the random seed to the current time.
 * Also takes in an argument for the number of points to use.
 * If not provided, the default number of points is 10000.
 * Reads the list of points from the default file.
 */
int main(int argc, const char* argv[]) {
  std::srand(time(NULL));

  if (argc > 1) {
    int npoints = atoi(argv[1]);
    std::vector<Point> points = generatePoints(npoints);
    savePoints(points);
  }

  auto points = readPoints();

  std::ofstream outfile("results.txt");
  timer(points, outfile, part2, "Recursive", 10);
  outfile.close();
}

#endif