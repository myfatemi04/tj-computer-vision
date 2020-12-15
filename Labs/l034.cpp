#ifndef LAB_PART
#define LAB_PART 4
#endif

#ifndef LAB_L034
#define LAB_L034

#include <iostream>
#include <vector>
#include <unordered_map>
#include "l03core.cpp"
#include "l033.cpp"

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

/**
 * Returns a closest point, at a pair, where the first point is the point found,
 * and the second point is the original point
 */
PointPair getClosestPoint(const GridMap& gridmap, const Point& point, double delta) {
	GridSquare center = makeGridSquare(point, delta);
	num startX = (center.first > 2) * (center.first - 2);
	num startY = (center.second > 2) * (center.second - 2);
	PointPair pair;

	for (num x = startX; x < startX + 5; x++) {
		for (num y = startY; y < startY + 5; y++) {
			GridSquare square(x, y);
			if (gridmap.count(square)) {
				pair.minify(PointPair(gridmap.at(square), point));
			}
		}
	}
	
	return pair;
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
	knuthShuffle(points);

	GridMap grid;
	PointPair closest(points[0], points[1]);
	
	for (size_t index = 0; index < points.size(); index++) {
		PointPair closestPoint = getClosestPoint(grid, points[index], closest.getDistance());

		if (closest.minify(closestPoint)) {
			grid.clear();

			for (size_t i = 0; i < index; i++) {
				grid[makeGridSquare(points[i], closest.getDistance())] = points[i];
			}
		}

		grid[makeGridSquare(points[index], closest.getDistance())] = points[index];
	}

	return closest;
}

#endif

#if LAB_PART == 4

void timeMany();

int main() {
	std::srand(time(NULL));
	// timeMany();

  auto points = readPoints("points100k.txt");
  // auto pointsDuplicate = std::vector<Point>(points);

  std::ofstream outfile("results.txt");
  // timer(pointsDuplicate, outfile, part3, "Recursive Optimized", 2);
  timer(points, outfile, part4, "Hashing Method", 2);
  outfile.close();
}

void timeMany() {
  std::ofstream outfile("results.txt");
  for (int i = 0; i < 12; i++) {
		auto points = readPoints(files[i]);
    timer(points, outfile, part4, "Dictionary", 10);
  }
	outfile.close();
}

#endif
