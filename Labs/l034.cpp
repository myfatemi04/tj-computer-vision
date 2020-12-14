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

num max(num a, num b) {
	return a > b ? a : b;
}

Point* getClosestPoint(const GridMap& gridmap, const Point& point, double delta) {
	auto center = makeGridSquare(point, delta);
	auto startX = (center.first > 2) * (center.first - 2);
	auto startY = (center.second > 2) * (center.second - 2);

	num closestX = 0;
	num closestY = 0;
	double closestDistance = 2;
	bool found = false;
	static GridSquare square(0, 0);

	for (auto x = startX; x < startX + 5; x++) {
		for (auto y = startY; y < startY + 5; y++) {
			square.first = x;
			square.second = y;
			if (gridmap.count(square)) {
				// If they are closer than 'delta'
				double distance = point.distance(gridmap.at(square));
				if (distance < delta) {
					if (!found) {
						// Set the initial value for 'closest'
						closestX = x;
						closestY = y;
						found = true;
					} else if (distance < closestDistance) {
						// Set the updated value for 'closest'
						closestX = x;
						closestY = y;
						closestDistance = distance;
					}
				}
			}
		}
	}
	
	if (found) {
		return new Point(gridmap.at(GridSquare(closestX, closestY)));
	}
	
	return nullptr;
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
	// std::cout << "shuffle:" << getMillis() - start << "ms\n";

	GridMap grid;
	Point *closestPoint1 = &points.at(0);
	Point *closestPoint2 = &points.at(1);
	double delta = closestPoint1->distance(*closestPoint2);

	long long timeSpentClearing = 0;
	long long loopStartTime = getMillis();
	for (int index = 0; index < points.size(); index++) {
		auto closestPoint = getClosestPoint(grid, points.at(index), delta);

		if (closestPoint != nullptr) {
			closestPoint1 = &points.at(index);
			closestPoint2 = closestPoint;

			delta = closestPoint1->distance(*closestPoint2);

			long long timeMarker1 = getMillis();
			grid.clear();
			
			long long timeMarker2 = getMillis();

			for (int i = 0; i < index; i++) {
				grid[makeGridSquare(points.at(i), delta)] = points.at(i);
			}
			
			long long timeMarker3 = getMillis();

			// std::cout << "[restart] grid clear " << (timeMarker2 - timeMarker1) << "ms, ";
			// std::cout << "grid remake " << (timeMarker3 - timeMarker2) << "ms\n";
			timeSpentClearing += (timeMarker3 - timeMarker1);
		}

		grid[makeGridSquare(points.at(index), delta)] = points.at(index);
	}

	long long loopEndTime = getMillis();
	long long timeSpentQuerying = loopEndTime - loopStartTime - timeSpentClearing;

	// std::cout << "Spent a total of " << timeSpentClearing << "ms remaking the grid\n";
	// std::cout << "Spent a total of " << timeSpentQuerying << "ms finding points (";
	// std::cout << (timeSpentQuerying / (double)points.size()) << "ms / point)\n";

	return PointPair(*closestPoint1, *closestPoint2);
}

#endif

#if LAB_PART == 4

int main() {
	std::srand(time(NULL));

  auto points = readPoints("points100k.txt");
  auto pointsDuplicate = std::vector<Point>(points);

  std::ofstream outfile("results.txt");
  timer(pointsDuplicate, outfile, part3, "Recursive Optimized", 2);
  timer(points, outfile, part4, "Hashing Method", 2);
  outfile.close();
}

#endif
