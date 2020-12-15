#ifndef LAB_PART
#define LAB_PART 3
#endif

#ifndef LAB_L033
#define LAB_L033

#include <algorithm>
#include <vector>
#include "l03core.cpp"
// #include "l032.cpp"

/**
 * Improved divide and conquer method. Optimizes based on the fact
 *  that you can make the comparisons between strips work in O(n log n)
 *  time, which is an improvement over O(n) time.
 */
PointPair helper3(std::vector<Point>& points, int begin, int end) {
  int numPoints = (end - begin);
  int mid = (begin + end) >> 1;
  if (numPoints == 2) {
    return PointPair(points[begin], points[end - 1]);
  } else if (numPoints == 3) {
    PointPair res;
    res.minify(PointPair(points[begin + 0], points[begin + 1]));
    res.minify(PointPair(points[begin + 1], points[begin + 2]));
    res.minify(PointPair(points[begin + 2], points[begin + 0]));

    return res;
  } else {
    PointPair closest;
    // Right side
    closest.minify(helper3(points, begin, mid));
    // Left side
    closest.minify(helper3(points, mid, end));
    // Closest pair where both are in left, or both are in right
    double d = closest.getDistance();

    // Maybe one of the points is on the left side and one is on the right?
    // Make a 'strip' of points down the middle
    int stripLeft = mid, stripRight = mid;
    double middleX = points[mid].getX();

    // [begin, mid) = left side
    while ((stripLeft > begin) && (points[stripLeft - 1].getX() >= (middleX - d))) {
      stripLeft--;
    }

    // [mid, end] = right side
    while ((stripRight < end) && (points[stripRight].getX() <= (middleX + d))) {
      stripRight++;
    }

    int stripSize = stripRight - stripLeft;

    // no point in sorting if you're checking the next 16 points anyway
    // otherwise, do the other method, as we can avoid sorting it again
    if (stripSize > 16) {
      // Create a strip sorted by y values
      // Here, we make a vector of pointers to points in the original vector.
      // This is so we do not lose the reference, which happens if you create
      // the strip based on the pointers to the first and last elements
      std::vector<Point> strip;
      for (int i = stripLeft; i < stripRight; i++) {
        strip.push_back(points[i]);
      }
      // Sorts the points
      std::sort(strip.begin(), strip.end(), comparePointYValues);

      // For each point in the strip, compare it with the next points
      // until the difference in y is greater than d
      for (int i = 0; i < stripSize; i++) {
        double maxY = strip[i].getY() + d;
        for (int j = i + 1; j < stripSize && strip[j].getY() < maxY; j++) {
          closest.minify({strip[i], strip[j]});
        }
      }
    } else {
      // For each point in the left side of the strip, compare it with
      // a point on the right side of the strip until the difference in
      // X is greater than d.
      for (int i = stripLeft; i < mid; i++) {
        double maxX = points[i].getX() + d;
        for (int j = mid; j < stripRight && points[j].getX() < maxX; j++) {
          closest.minify({points[i], points[j]});
        }
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

#endif

#if LAB_PART == 3

void timeMany();

/**
 * Main method. Sets the random seed to the current time.
 * Also takes in an argument for the number of points to use.
 * If not provided, the default number of points is 10000.
 * Reads the list of points from the default file.
 */
int main(int argc, const char* argv[]) {
  std::srand(time(NULL));
  timeMany();

  // std::vector<Point> points;
  // if (argc > 1) {
  //   points = readPoints(argv[1]);
  // } else {
  //   points = readPoints();
  // }

  // std::ofstream outfile("results.txt");
  // timer(points, outfile, part3, "Recursive Optimized", 10);
}

void timeMany() {
  std::ofstream outfile("results.txt");
  for (int i = 0; i < 12; i++) {
    auto points = readPoints(files[i]);
    timer(points, outfile, part3, "Recursive Optimized", 10);
  }
  outfile.close();
}

#endif