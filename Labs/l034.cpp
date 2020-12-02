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
		bool addPoint(Point& point) {
			GridSquare gs = this -> getGridSquare(point);
		}

		Point getPointAt(GridSquare& gridSquare) const {
			return gridSquares.at(gridSquare);
		}
		
		bool hasPointAt(GridSquare& gridSquare) const {
			return (gridSquares.find(gridSquare) != gridSquares.end());
		}

		/**
		 * Looks for the closest point within the 5x5 grid of neighboring grid squares. If no points are found in those grid squares, this returns -1.
		 */
		PointPair* getClosestPointWithinBoundaries(Point& point) {
			PointPair *closest = NULL;
			GridSquare gridSquare = this -> getGridSquare(point);
			for (int x = -2; x <= 2; x++) {
				for (int y = -2; y <= 2; y++) {
					GridSquare gridSquareToCheck = GridSquare(x, y);
					if (this -> hasPointAt(gridSquareToCheck)) {
						Point pointToCheck = this -> getPointAt(gridSquareToCheck);
						double distance = Point::distance(pointToCheck, point);
						if (closest == NULL || distance < closest -> getDistance()) {
							if (closest) {
								delete closest;
							}
							closest = new PointPair(point, pointToCheck);
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
		Point p = points.at(i);
		PointPair* closest = grid.getClosestPointWithinBoundaries(p);
		if (closest != NULL) {
			if (closest -> getDistance() < grid.getDelta()) {
				grid.clear();
				grid.setDelta(closest -> getDistance());

				for (int j = 0; j <= i; j++) {
					Point pointToAdd = points.at(j);
					grid.addPoint(pointToAdd);
				}
			}
		}
	}
}

#endif
