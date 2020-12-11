#ifndef LAB_PART
#define LAB_PART 4
#endif

#ifndef LAB_L034
#define LAB_L034

#include <iostream>
#include <vector>
#include <unordered_map>
#include "l03core.cpp"

typedef unsigned long long num;
typedef std::pair<num, num> GridSquare;

struct hash_pair {
	template <class T1, class T2>
	num operator()(const std::pair<T1, T2>& p) const {
		auto hash1 = std::hash<T1>{}(p.first);
		auto hash2 = std::hash<T2>{}(p.second);
		return hash1 ^ hash2;
	}
};

typedef std::unordered_map<GridSquare, Point, hash_pair> GridMap;

GridSquare makeGridSquare(const Point& point, double delta) {
	return GridSquare(
		(num)(point.getX() / (delta / 2)),
		(num)(point.getY() / (delta / 2))
	);
}

num max(num a, num b) {
	return a > b ? a : b;
}

Point* getClosestPoint(const GridMap& gridmap, const Point& point, double delta) {
	auto square = makeGridSquare(point, delta);
	auto startX = (square.first > 2) * (square.first - 2);
	auto startY = (square.second > 2) * (square.second - 2);

	Point closest;
	bool found = false;

	for (auto x = startX; x < startX + 4; x++) {
		for (auto y = startY; y < startY + 4; y++) {
			if (gridmap.count(GridSquare(x, y)) > 0) {
				auto point_ = gridmap.at(GridSquare(x, y));
				// If they are closer than 'delta'
				double distance = point.distance(point_);
				if (distance < delta) {
					if (!found) {
						// Set the initial value for 'closest'
						closest = point_;
						found = true;
					} else if (distance < closest.distance(point)) {
						// Set the updated value for 'closest'
						closest = point_;
					}
				}
			}
		}
	}

	if (found) {
		return new Point(closest);
	} else {
		return nullptr;
	}
}

void knuthShuffle(std::vector<Point>& points) {
	for (int currentIndex = points.size() - 1; currentIndex >= 0; currentIndex--) {
		// each time, choose a point to move to the back of the vector
		// it must be before the back of the vector
		int swapIndex = (int)(getRandom() * currentIndex);

		// swap index i and the random index
		auto currentPoint = points[currentIndex];
		auto swapPoint = points[swapIndex];

		points[currentIndex] = swapPoint;
		points[swapIndex] = currentPoint;
	}
}

PointPair part4(std::vector<Point>& points) {
	long long start = getMillis();
	knuthShuffle(points);
	std::cout << "shuffle:" << getMillis() - start << "ms\n";

	GridMap grid;
	Point *closestPoint1 = &points.at(0);
	Point *closestPoint2 = &points.at(1);
	double delta = closestPoint1->distance(*closestPoint2);
	
	for (int index = 0; index < points.size(); index++) {
		auto closestPoint = getClosestPoint(grid, points.at(index), delta);

		if (closestPoint != nullptr) {
			closestPoint1 = &points.at(index);
			closestPoint2 = closestPoint;

			delta = closestPoint1->distance(*closestPoint2);
			grid.clear();
			for (int i = 0; i < index; i++) {
				grid[makeGridSquare(points.at(i), delta)] = points.at(i);
			}
		}

		grid[makeGridSquare(points[index], delta)] = points[index];
	}

	return PointPair(*closestPoint1, *closestPoint2);
}

#endif

#if LAB_PART == 4

int main() {
	std::srand(time(NULL));

  auto points = readPoints("points1m.txt");

  std::ofstream outfile("results.txt");
  timer(points, outfile, part4, "Hashing Method", 2);
  outfile.close();
}

#endif
