#ifndef LAB_PART
#define LAB_PART 1
#endif

#ifndef LAB_L031
#define LAB_L031

#include <iomanip>
#include <vector>
#include "l03core.cpp"

PointPair part1(std::vector<Point>& points) {
  PointPair closest(points.at(0), points.at(1));

  for (auto i = points.begin(); i != points.end(); i++) {
    for (auto j = std::next(i); j != points.end(); j++)  {
      closest.minify({*i, *j});
    }
  }

  return closest;
}

#endif

#if LAB_PART == 1
int main(int argc, const char* argv[]) {
  std::srand(time(NULL));

  if (argc > 1) {
    std::vector<Point> points = generatePoints(atoi(argv[1]));
    savePoints(points);
  }

  auto points = readPoints();

  std::ofstream outfile("results.txt");
  timer(points, outfile, part1, "Brute Force", 10);
  outfile.close();
}

#endif