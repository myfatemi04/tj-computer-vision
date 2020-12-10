#ifndef LAB_PART
#define LAB_PART 4
#endif

#ifndef LAB_L034
#define LAB_L034

#include <vector>
#include <unordered_map>
#include "geometry.cpp"
#include "l03core.cpp"

typedef unsigned long long num;
typedef std::pair<num, num> GridSquare;

struct hash_pair {
	template <class T1, class T2>
	size_t operator()(const std::pair<T1, T2>& p) const {
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

GridMap* makeGridMap(
	std::vector<Point>::iterator start,
	std::vector<Point>::iterator end,
	double delta
) {
	GridMap *gridMap = new GridMap();

	for (auto i = start; i != end; i++) {
		if (gridMap -> count(makeGridSquare(*i, delta))) {
			std::cerr << delta << "Already exists- " << *i << '\n';
			exit(0);
		}

		gridMap -> emplace(makeGridSquare(*i, delta), *i);
	}

	return gridMap;
}

num max(num a, num b) {
	return a > b ? a : b;
}

Point* getClosestPoint(const GridMap& gridmap, const Point& point, double delta) {
	auto square = makeGridSquare(point, delta);
	auto startX = square.first > 2 ? square.first - 2 : 0;
	auto startY = square.second > 2 ? square.second - 2 : 0;

	Point* closest = nullptr;

	for (auto x = startX; x < startX + 4; x++) {
		for (auto y = startY; y < startY + 4; y++) {
			if (gridmap.count(GridSquare(x, y)) > 0) {
				auto point_ = gridmap.at(GridSquare(x, y));
				// If they are closer than 'delta'
				double distance = Point::distance(point_, point);
				if (distance < delta) {
					if (closest == nullptr) {
						// Set the initial value for 'closest'
						closest = new Point(point_);
					} else if (distance < Point::distance(*closest, point)) {
						// Set the updated value for 'closest'
						delete closest;
						closest = new Point(point_);
					}
				}
			}
		}
	}

	return closest;
}

class Grid {
	private:
		GridMap* gridSquares;
		double delta = 1;
		Point *closest1, *closest2;

		void updateClosestPair(Point* closest1, Point *closest2) {
			this -> closest1 = closest1;
			this -> closest2 = closest2;
			// Essentially, recreate the grid with a new 'delta' value
			// 'delta' value is the width of each grid square
			this -> setDelta(Point::distance(*closest1, *closest2));
		}

	public:
		Grid() {
			this -> delta = 1;
			this -> gridSquares = new GridMap();
		}

		num getGridX(const Point& p) const noexcept {
			return (num)(p.getX() / (this -> delta / 2));
		}

		num getGridY(const Point& p) const noexcept {
			return (num)(p.getY() / (this -> delta / 2));
		}

		GridSquare getGridSquare(const Point& p) const noexcept {
			return GridSquare(getGridX(p), getGridY(p));
		}

		/**
		 * Adds a point to the grid at the correct grid area.
		 * If the point already exists returns false.
		 */
		void addPoint(const Point& point) {
			Point* pointWithinDelta = this -> getPointWithinDelta(point);
			if (pointWithinDelta != nullptr) {
				double distance = Point::distance(point, *pointWithinDelta);
				this -> updateClosestPair(new Point(point), pointWithinDelta);
			}
			delete pointWithinDelta;

			gridSquares -> emplace(getGridSquare(point), point);
		}

		Point getPointForGridSquare(const GridSquare& gridSquare) const {
			return gridSquares -> at(gridSquare);
		}
		
		bool hasPointAt(const GridSquare& gridSquare) const {
			return gridSquares -> count(gridSquare);
		}

		/**
		 * Looks for the closest point within the 5x5 grid of neighboring grid squares. If no points are found in those grid squares, this returns -1.
		 */
		Point* getPointWithinDelta(const Point& point) const {
			Point* closest = nullptr;

			num centerX = getGridX(point), centerY = getGridY(point);
			num startX = 0, startY = 0;
			if (centerX > 2) startX = centerX - 2;
			if (centerY > 2) startY = centerY - 2;

			// Iterate over nearest grid squares
			for (num x = startX; x <= centerX + 2; x++) {
				for (num y = startY; y <= centerY + 2; y++) {
					GridSquare gridSquare(x, y);
					if (!hasPointAt(gridSquare)) {
						continue;
					}

					Point pointToCheck = this -> gridSquares -> at(gridSquare);
					double distance = Point::distance(point, pointToCheck);
					if (distance >= this -> delta) {
						continue;
					}

					if (closest == nullptr) {
						closest = new Point(pointToCheck);
					} else if (distance < Point::distance(point, *closest)) {
						delete closest;
						closest = new Point(pointToCheck);
					}
				}
			}

			return closest;
		}

		int size() const noexcept {
			return this -> gridSquares -> size();
		}

		void setDelta(double newDelta) {
			GridMap* oldGridSquares = this -> gridSquares;
			this -> gridSquares = new GridMap();
			this -> delta = newDelta;
			for (auto pair: *oldGridSquares) {
				this -> addPoint(pair.second);
			}
			delete oldGridSquares;
		}

		double getDelta() const {
			return this -> delta;
		}

		PointPair* getClosestPair() const {
			return new PointPair(*this -> closest1, *this -> closest2);
		}
};

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

	GridMap* grid = new GridMap();
	Point *closestPoint1 = &points.at(0);
	Point *closestPoint2 = &points.at(1);
	double delta = 1;
	
	for (auto it = points.begin(); it != points.end(); it = std::next(it)) {
		auto closestPoint = getClosestPoint(*grid, *it, delta);

		if (closestPoint != nullptr) {
			closestPoint1 = new Point(*it);
			closestPoint2 = new Point(*closestPoint);
			delta = Point::distance(*closestPoint, *it);

			delete grid;
			grid = makeGridMap(points.begin(), it, delta);
		}

		grid -> emplace(makeGridSquare(*it, delta), *it);
	}

	return PointPair(*closestPoint1, *closestPoint2);
}

#endif

#if LAB_PART == 4

int main() {
	std::srand(time(NULL));

  auto points = readPoints();

  std::ofstream outfile("results.txt");
  timer(points, outfile, part4, "Hashing Method", 1);
  outfile.close();
}

#endif
