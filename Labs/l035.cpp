#include "l03core.cpp"
#include <vector>
#include <iostream>
#include <iomanip>

const int MAX_POINTS = 16;

double getMedian(const double *array, const int length) {
	int end = length - 1;
	// length / 2 = floor(length / 2)
	int lowerMiddle = end / 2;
	// length / 2 + length % 2 = ceil(length / 2)
	int upperMiddle = lowerMiddle + end & 1;

	return (array[lowerMiddle] + array[upperMiddle]) / 2.0;
}

bool between(double a, double b, double x) {
	return (a < x) && (x < b);
}

class Rectangle {
	private:
		double minX, maxX, minY, maxY;

	public:
		Rectangle(
			double minX,
			double maxX,
			double minY,
			double maxY
		):
			minX(minX),
			maxX(maxX),
			minY(minY),
			maxY(maxY) {}

		bool fullyContains(const Rectangle& other) const {
			return (
				(
					between(minX, maxX, other.minX) &&
					between(minX, maxX, other.maxX)
				) &&
				(
					between(minY, maxY, other.minY) &&
					between(minY, maxY, other.maxY)
				)
			);
		}

		bool doesNotIntersect(const Rectangle& other) const {
			return (
				(minX > other.maxX || maxX < other.minX) ||
				(minY > other.maxY || maxY < other.minY)
			);
		}

		bool intersects(const Rectangle& other) {
			return !doesNotIntersect(other);
		}

		bool contains(const Point& point) const {
			return (
				between(minX, maxX, point.getX()) &&
				between(minY, maxY, point.getY())
			);
		}

		double getMinX() const {
			return minX;
		}

		double getMinY() const {
			return minY;
		}

		double getMaxX() const {
			return maxX;
		}

		double getMaxY() const {
			return maxY;
		}

		static Rectangle getRandomRectangle() {
			double a = getRandom();
			double b = getRandom();
			double c = getRandom();
			double d = getRandom();
			return Rectangle(a, a + b, c, c + d);
		}
};

std::ostream& operator<<(std::ostream& stream, const Rectangle& out) {
	return 
		stream
			<< out.getMinX() << ">" << out.getMaxX() << ", "
			<< out.getMinY() << ">" << out.getMaxY();
}

PointPair closestPairOfPointsBruteForce(const std::vector<Point>& points) {
	// Brute force search
	PointPair closest = PointPair(points[0], points[1]);

	for (int i = 0; i < points.size(); i++) {
		for (int j = i + 1; j < points.size(); j++) {
			closest.minify(PointPair(points[i], points[j]));
		}
	}
	
	return closest;
}

PointPair closestPairOfPointsStrip(std::vector<Point>& points) {
	PointPair closestPair;
	for (int i = 0; i < points.size(); i++) {
		int maxJ = i + 8 < points.size() ? i + 8 : points.size();
		for (int j = i; j < maxJ; j++) {
			closestPair.minify(PointPair(points.at(i), points.at(j)));
		}
	}
	return closestPair;
}

const int X_BOTTOM_HALF = 0, Y_BOTTOM_HALF = 0, X_TOP_HALF = 1, Y_TOP_HALF = 2;

class TreeNode {
	private:
		double midX, midY;
		Rectangle boundaries;
		bool divided = false;

		// list of the branches
		// two bits: the first bit is if point > y value,
		// the second bit is if point > x value.
		// example: 0b10 [2] means point > y critical value but point <= x critical value
		// example: 0b00 [0] means point < both critical values
		// example: 0b11 [3] means point > both critical values
		TreeNode *branches;
		std::vector<Point> points;

		/**
		 * The "Optimal Split" is the one that creates the greatest separation
		 * Additionally, by creating separation, it splits based on which direction
		 * has the higher distance between the highest and lowest points
		 */
		void createOptimalSplit() {
			double xMax, xMin, yMax, yMin;
			xMin = xMax = points[0].getX();
			yMin = yMax = points[0].getY();

			for (Point p : points) {
				xMin = p.getX() < xMin ? p.getX() : xMin;
				xMax = p.getX() > xMax ? p.getX() : xMax;
				yMin = p.getY() < yMin ? p.getY() : yMin;
				yMax = p.getY() > yMax ? p.getY() : yMax;
			}

			divided = true;
			midX = (xMin + xMax) / 2;
			midY = (yMin + yMax) / 2;
			// midX = (boundaries.getMinX() + boundaries.getMaxX()) / 2;
			// midY = (boundaries.getMinY() + boundaries.getMaxY()) / 2;

			// Create the subtrees
			branches = new TreeNode[4] {
				// low low
				TreeNode(
					Rectangle(
						boundaries.getMinX(),
						midX,
						boundaries.getMinY(),
						midY
					)
				),
				// high low
				TreeNode(
					Rectangle(
						midX,
						boundaries.getMaxX(),
						boundaries.getMinY(),
						midY
					)
				),
				// low high
				TreeNode(
					Rectangle(
						boundaries.getMinX(),
						midX,
						midY,
						boundaries.getMaxY()
					)
				),
				// high high
				TreeNode(
					Rectangle(
						midX,
						boundaries.getMaxX(),
						midY,
						boundaries.getMaxY()
					)
				),
			};

			// Re-add the points
			for (Point p : points) {
				addPointToCorrectSubtree(p);
			}
		}

		void addPointToCorrectSubtree(const Point& point) {
			branches[getSubtree(point)].addPoint(point);
		}
	
	public:
		TreeNode(Rectangle boundaries): boundaries(boundaries) {}

		int getSubtree(const Point& point) {
			return
				((point.getX() > midX) << 0) |
				((point.getY() > midY) << 1);
		}

		/**
		 * Adds a point to this TreeNode. If the TreeNode is at
		 * max capacity, it splits into to other TreeNodes, finding the optimal
		 * split to do so.
		 */
		void addPoint(const Point& point) {
			if (!divided) {
				if (points.size() == MAX_POINTS) {
					createOptimalSplit();
				} else {
					points.push_back(point);
					return;
				}
			}
			
			addPointToCorrectSubtree(point);
		}
		
		/**
		 * Returns the value used when choosing a subtree to place the next point into.
		 * Will be an unassigned value if the node has not been divided yet.
		 * You can check if it has been divided or not with TreeNode::getDivideType()
		 */
		double getXSplit() const {
			return midX;
		}

		double getYSplit() const {
			return midY;
		}

		bool isDivided() const {
			return divided;
		}

		int getPointCount() const {
			if (!divided) {
				return points.size();
			} else {
				int count = 0;
				for (int i = 0; i < 4; i++) {
					count += branches[i].getPointCount();
				}
				return count;
			}
		}

		/**
		 * Finds all points within a certain radius of another point.
		 * Radius inclusive.
		 */
		void addPointsWithinRadius(
			std::vector<Point>& out,
			const Point& center,
			double radius
		) const {
			// If we don't intersect the circumscribed rectangle,
			// we are guaranteed not to have any points within the circle
			Rectangle circumscribedRectangle(
				center.getX() - radius,
				center.getX() + radius,
				center.getY() - radius,
				center.getY() + radius
			);

			const double sqrt2 = 1.41421356;

			Rectangle inscribedRectangle(
				center.getX() - radius * sqrt2 / 2,
				center.getX() + radius * sqrt2 / 2,
				center.getY() - radius * sqrt2 / 2,
				center.getY() + radius * sqrt2 / 2
			);

			if (circumscribedRectangle.intersects(boundaries)) {
				if (divided) {
					for (int i = 0; i < 4; i++) {
						branches[i].addPointsWithinRadius(out, center, radius);
					}
				} else if (inscribedRectangle.fullyContains(boundaries)) {
					out.insert(out.end(), points.begin(), points.end());
				} else {
					for (Point point : points) {
						if (Point::areWithinDistance(center, point, radius)) {
							out.push_back(point);
						}
					}
				}
			}
		}

		Point* findClosestPoint(const Point& center, double radius) const {
			if (!divided) {
				Point closest;
				bool found = false;
				double closestDistance;

				for (Point p : points) {
					double distance = Point::distance(center, p);
					if (distance <= radius) {
						if (!found) {
							found = true;
							closest = p;
							closestDistance = distance;
						} else if (distance < closestDistance) {
							closest = p;
							closestDistance = distance;
						}
					}
				}

				if (found) {
					return new Point(closest);
				}
			} else {
				// If we don't intersect the circumscribed rectangle,
				// we are guaranteed not to have any points within the circle
				Rectangle circumscribedRectangle(
					center.getX() - radius,
					center.getX() + radius,
					center.getY() - radius,
					center.getY() + radius
				);

				if (circumscribedRectangle.intersects(boundaries)) {
					Point closest;
					bool found = false;
					double closestDistance;

					for (int i = 0; i < 4; i++) {
						Point *closestThisBranch = branches[i].findClosestPoint(center, radius);
						if (closestThisBranch != nullptr) {
							double distance = Point::distance(center, *closestThisBranch);
							if (found) {
								if (distance < closestDistance) {
									closestDistance = distance;
									closest = *closestThisBranch;
								}
							} else {
								found = true;
								closest = *closestThisBranch;
								closestDistance = distance;
							}
							delete closestThisBranch;
						}
					}

					if (found) {
						return new Point(closest);
					}
				}
			}

			return nullptr;			
		}

		/**
		 * Finds all points within a certain X and Y boundary
		 */
		void addPointsInRectangle(
			std::vector<Point>& out,
			const Rectangle& rectangle
		) const {
			if (rectangle.doesNotIntersect(boundaries)) {
				return;
			}

			if (!divided) {
				// If fully contained by rect, add all points
				if (rectangle.fullyContains(boundaries)) {
					out.insert(out.end(), points.begin(), points.end());
				} else {
					for (Point p : points) {
						if (rectangle.contains(p)) {
							out.push_back(p);
						}
					}
				}
			} else {
				for (int i = 0; i < 4; i++) {
					branches[i].addPointsInRectangle(out, rectangle);
				}
			}
		}

		PointPair findClosestPairOfPoints() const {
			// Run different code based on if this tree is split or not
			if (!divided) {
				return closestPairOfPointsBruteForce(points);
			} else {
				// Find the closest pairs in each subtree
				PointPair closestPair;
				for (int i = 0; i < 4; i++) {
					if (branches[i].getPointCount() >= 2) {
						closestPair.minify(branches[i].findClosestPairOfPoints());
					}
				}

				// Find the closest pairs in bordering subtrees
				std::vector<Point> strip;

				// Vertical strip
				addPointsInRectangle(
					strip,
					Rectangle(
						midX - closestPair.getDistance(),
						midX + closestPair.getDistance(),
						boundaries.getMinY(),
						boundaries.getMaxY()
					)
				);
				if (strip.size() >= 2) {
					// std::sort(strip.begin(), strip.end(), comparePointYValues);
					closestPair.minify(closestPairOfPointsBruteForce(strip));
				}

				// Horizontal strip
				strip.clear();
				addPointsInRectangle(
					strip,
					Rectangle(
						boundaries.getMinX(),
						boundaries.getMaxX(),
						midY - closestPair.getDistance(),
						midY + closestPair.getDistance()
					)
				);
				if (strip.size() >= 2) {
					std::sort(strip.begin(), strip.end(), comparePointXValues);
					closestPair.minify(closestPairOfPointsBruteForce(strip));
				}

				return closestPair;
			}
		}

		void addAllPoints(std::vector<Point>& points) const {
			if (divided) {
				for (int i = 0; i < 4; i++) {
					branches[i].addAllPoints(points);
				}
			} else {
				points.insert(points.end(), this->points.begin(), this->points.end());
			}
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

int main() {
	std::srand(time(NULL));
	std::cout << std::setprecision(17);

	std::cout << "Started program\n";

	std::vector<Point> points = readPoints("points100k_random.txt");
	
	std::cout << "Read points\n";

	knuthShuffle(points);

	std::cout << "Shuffled points\n";

	long long startTime = getMillis();

	TreeNode tree(Rectangle(0, 1, 0, 1));

	/*
	PointPair closestPair(points[0], points[1]);

	for (Point point : points) {
		Point *closest = tree.findClosestPoint(point, closestPair.getDistance());
		if (closest != nullptr) {
			closestPair.minify(PointPair(point, *closest));
		}
		tree.addPoint(point);
	}

	std::cout << closestPair;
	*/

	for (Point point : points) {
		tree.addPoint(point);
	}

	long long endTime = getMillis();
	
	std::cout << "Created tree in " << (endTime - startTime) << "ms\n";

	startTime = getMillis();
	PointPair closest;
	int cycles = 10;
	for (int i = 0; i < cycles; i++) {
		closest = tree.findClosestPairOfPoints();
	}
	endTime = getMillis();

	std::cout << "Optimized 2D tree x [" << cycles << "] " << (endTime - startTime) << "ms ";
	std::cout << "(" << (endTime - startTime)/ cycles << "/cycle)\n";
	std::cout << closest << '\n';
}

void benchmarkRectangles(const std::vector<Point>& points, const TreeNode& tree) {
	std::vector<Point> out;

	long long startTime, endTime;
	double x, y, r;

	startTime = getMillis();

	for (int i = 0; i < 500; i++) {
		// auto rect = Rectangle::getRandomRectangle();
		x = getRandom();
		y = getRandom();
		r = getRandom() / 2;
		out.clear();
		tree.addPointsWithinRadius(out, Point(x, y), r);
		// tree.addPointsInRectangle(out, rect);
	}
	endTime = getMillis();

	std::cout << "Found in " << (endTime - startTime) << "ms w/ tree\n";

	startTime = getMillis();
	for (int i = 0; i < 500; i++) {
		// auto rect = Rectangle::getRandomRectangle();
		x = getRandom();
		y = getRandom();
		r = getRandom() / 2;
		out.clear();
		for (Point point : points) {
			// if (rect.contains(point)) {
			if (Point::areWithinDistance(point, Point(x, y), r)) {
				out.push_back(point);
			}
		}
	}
	endTime = getMillis();

	std::cout << "Found in " << (endTime - startTime) << "ms w/ brute force\n";
}

// Now, we can use this data to our advantage
// Find the closest pair of points in a subtree
