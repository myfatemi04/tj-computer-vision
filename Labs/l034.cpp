#ifndef LAB_PART
#define LAB_PART 4
#endif

#ifndef LAB_L034
#define LAB_L034

#include <vector>
#include <unordered_map>
#include "geometry.cpp"
#include "l03core.cpp"

class GridSquare {
	private:
		int x, y;

	public:
		GridSquare(int x, int y): x(x), y(y) {}

		int getX() const { return x; }
		int getY() const { return y; }

		bool operator==(const GridSquare& other) const {
			return other.getX() == this->getX() && other.getY() == this->getY();
		}
};

namespace std {
	template<>
	struct hash<GridSquare> {
		size_t operator()(const GridSquare& gridSquare) const {
			hash<int> intHash = hash<int>();
			return intHash(gridSquare.getX()) ^ (intHash(gridSquare.getY()) << 1);
		}
	};
}

class Grid {
	private:
		std::unordered_map<GridSquare, Point> gridSquares;
		double delta;

	public:
		Grid(double initialDelta): delta(initialDelta) {}
		~Grid() {
			delete &gridSquares;
		}

		GridSquare getGridSquare(Point p) const {
			int x = (int)(p.getX() / delta);
			int y = (int)(p.getY() / delta);
			return GridSquare(x, y);
		}

		/**
		 * Adds a point to the grid at the correct grid area.
		 * If the point already exists returns false.
		 */
		void addPoint(const Point& point) {
			gridSquares.emplace(getGridSquare(point), point);
		}

		Point* getPointAt(GridSquare& gridSquare) const {
			if (this -> hasPointAt(gridSquare)) {
				return new Point(gridSquares.at(gridSquare));
			} else {
				return NULL;
			}
		}
		
		bool hasPointAt(GridSquare& gridSquare) const {
			return (gridSquares.find(gridSquare) != gridSquares.end());
		}

		/**
		 * Looks for the closest point within the 5x5 grid of neighboring grid squares. If no points are found in those grid squares, this returns -1.
		 */
		PointPair* getClosestPointWithinBoundaries(Point& point) const {
			PointPair *closest = NULL;
			GridSquare gridSquare = this -> getGridSquare(point);

			// Iterate over nearest grid squares
			int centerX = gridSquare.getX();
			int centerY = gridSquare.getY();

			for (int x = centerX - 2; x <= centerX + 2; x++) {
				for (int y = centerY - 2; y <= centerY + 2; y++) {
					GridSquare gridSquareToCheck = GridSquare(x, y);
					Point *pointToCheck = this -> getPointAt(gridSquareToCheck);

					if (pointToCheck != NULL) {
						double distance = Point::distance(*pointToCheck, point);

						if (closest == NULL) {
							closest = new PointPair(point, *pointToCheck);
						} else {
							if (distance <closest -> getDistance()) {
								delete closest;
								closest = new PointPair(point, *pointToCheck);
							}
						}
					}
				}
			}

			return closest;
		}

		void clear() {
			gridSquares.clear();
		}

		void setDelta(double newDelta) {
			this -> delta = newDelta;
		}

		double getDelta() {
			return this -> delta;
		}
};

PointPair part4(std::vector<Point>& points) {
	PointPair closest = PointPair(points.at(0), points.at(1));
	Grid grid(closest.getDistance());
	
	for (int i = 0; i < points.size(); i++) {
		Point point = points.at(i);
		PointPair* closestHere = grid.getClosestPointWithinBoundaries(point);

		if (closestHere != NULL) {
			if (closestHere -> getDistance() < grid.getDelta()) {
				closest = *closestHere;

				grid.clear();
				grid.setDelta(closestHere -> getDistance());

				for (int j = 0; j <= i; j++) {
					grid.addPoint(points.at(j));
				}

				continue;
			}
		}

		grid.addPoint(point);
	}

	return closest;
}

#endif

#if LAB_PART == 4

int main() {
	std::srand(time(NULL));

  // if (argc > 1) {
  //   int npoints = atoi(argv[1]);
  //   std::vector<Point> points = generatePoints(npoints);
  //   savePoints(points);
  // }

  auto points = readPoints();

  std::ofstream outfile("results.txt");
  timer(points, outfile, part4, "Hashing Method", 10);
  outfile.close();
}

#endif
