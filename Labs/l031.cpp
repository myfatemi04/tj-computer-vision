#include<iomanip>
#include<list>
#include "geometry.cpp"
#include "image.cpp"
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

#ifndef LAB_03_MAIN
#define LAB_03_MAIN

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