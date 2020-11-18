#include<iomanip>
#include<list>
#include "geometry.cpp"
#include "image.cpp"
#include "l03core.cpp"

std::list<Point> readPointsList(int npoints) {
  std::ifstream in("points.txt");
  std::list<Point> points;
  for (int i = 0; i < npoints; i++) {
    double x, y;
    in >> x >> y;
    points.push_back(Point(x, y));
  }
  return points;
}

std::list<Point> generatePointsList(int npoints) {
  std::list<Point> points;
  for (int i = 0; i < npoints; i++) {
    points.push_back(Point(getRandom(), getRandom()));
  }
  return points;
}

void savePointsList(std::list<Point> points) {
  std::ofstream out("points.txt");
  
  // makes all numbers write with exactly 17 digits after the decimal point
  out << std::fixed << std::setprecision(17);
  for (auto it = points.begin(); it != points.end(); it++) {
    out << it->getX() << "  " << it->getY() << "\n";
  }
  
  out.close();
}

void part1() {
  std::list<Point> points = generatePointsList(50);
  std::ofstream outImage("points.ppm");
  savePointsList(points);

  // Closest points
  PointPair closest(*points.begin(), *(std::next(points.begin())));
  Image image(800, 800, 3, new int[3] {1, 1, 1});

  for (auto it1 = points.begin(); it1 != points.end(); it1++) {
    image.drawPointCircle(*it1 * 800, 2, new int[3] {0, 0, 0});
    for (auto it2 = std::next(it1); it2 != points.end(); it2++) {
      closest.minify({*it1, *it2});
    }
  }

  image.drawPointCircle(closest.getA() * 800, 2, new int[3] {1, 0, 0});
  image.drawPointCircle(closest.getB() * 800, 2, new int[3] {1, 0, 0});
  image.saveImage(outImage);
  outImage.close();
}

int main() {
  part1();
}