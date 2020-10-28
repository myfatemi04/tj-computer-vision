#include<iomanip>
#include<list>
#include "geometry.cpp"
#include "image.cpp"

std::list<Point> readPoints(std::istream& in, int npoints) {
  std::list<Point> points;
  for (int i = 0; i < npoints; i++) {
    double x, y;
    in >> x >> y;
    points.push_back(Point(x, y));
  }
  return points;
}

std::list<Point> generatePoints(int npoints) {
  std::list<Point> points;
  for (int i = 0; i < npoints; i++) {
    points.push_back(Point(getRandom(), getRandom()));
  }
  return points;
}

void savePoints(std::ostream& out, std::list<Point> points) {
  // makes all numbers write with exactly 17 digits after the decimal point
  out << std::fixed << std::setprecision(17);
  for (auto it = points.begin(); it != points.end(); it++) {
    out << it->getX() << "  " << it->getY() << "\n";
  }
}

int main() {
  std::list<Point> points = generatePoints(50);
  std::ofstream outPoints("points.txt");
  std::ofstream outImage("points.ppm");
  savePoints(outPoints, points);
  outPoints.close();

  // Closest points
  Point *p1, *p2;
  double minDistance = 2;
  Image image(800, 800, 3, new int[3] {1, 1, 1});

  for (auto it1 = points.begin(); it1 != points.end(); it1++) {
    image.drawPointCircle(*it1 * 800, new int[3] {0, 0, 0});
    for (auto it2 = std::next(it1); it2 != points.end(); it2++) {
      double dist = Point::distance(*it1, *it2);
      if (dist < minDistance) {
        minDistance = dist;
        p1 = &(*it1);
        p2 = &(*it2);
      }
    }
  }

  image.drawPointCircle(*p1 * 800, new int[3] {1, 0, 0});
  image.drawPointCircle(*p2 * 800, new int[3] {1, 0, 0});
  image.saveImage(outImage);
  outImage.close();
}