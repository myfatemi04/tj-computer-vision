#include <vector>
#include <iostream>

size_t max(size_t a, size_t b) {
	return (a > b) ? (a) : (b);
}

class Point;
class PPM;

typedef std::vector<Point> Points;
typedef Points Means;
typedef Points Cluster;
typedef std::vector<Cluster> Clusters;
typedef struct {
	int r, g, b;
} Pixel;

const Pixel black  = { 0x00, 0x00, 0x00 }; // 000
const Pixel blue   = { 0x00, 0x00, 0xFF }; // 001
const Pixel green  = { 0x00, 0xFF, 0x00 }; // 010
const Pixel cyan   = { 0x00, 0xFF, 0xFF }; // 011
const Pixel red    = { 0xFF, 0x00, 0x00 }; // 100
const Pixel purple = { 0xFF, 0x00, 0xFF }; // 101
const Pixel yellow = { 0xFF, 0xFF, 0x00 }; // 110
const Pixel white  = { 0xFF, 0xFF, 0xFF }; // 111

const Pixel orange = { 0xFF, 0x7F, 0x00 };

class Point {
	private:
		size_t dimensionCount;
		double *dimensions;

	public:
		Point(): dimensionCount(0), dimensions(nullptr) {}

		Point(double x, double y): dimensionCount(2), dimensions(new double[2] { x, y }) {}

		Point(double x, double y, double z): dimensionCount(3), dimensions(new double[3] { x, y, z }) {}

		Point(double *dimensions, int dimensionCount): dimensions(dimensions), dimensionCount(dimensionCount) {}

		int getDimensionCount() const {
			return dimensionCount;
		}

		double calculateDistanceSquared(const Point& other) const {
			double sum = 0;
			for (size_t dimension = 0; dimension < max(dimensionCount, other.dimensionCount); dimension++) {
				double difference = other[dimension] - at(dimension);
				sum += difference * difference;
			}
			return sum;
		}

		double at(size_t dimension) const {
			if (dimension < dimensionCount) {
				return dimensions[dimension];
			} else {
				return 0;
			}
		}

		double operator[](size_t dimension) const {
			return at(dimension);
		}

		bool operator==(const Point& other) const {
			if (dimensionCount != other.dimensionCount) {
				return false;
			}
			for (size_t i = 0; i < dimensionCount; i++) {
				if (at(i) != other.at(i)) {
					return false;
				}
			}
			return true;
		}

		/**
		 * Returns the index within the input vector of the closest mean
		 */
		int classify(const Means& means) const {
			int closestIndex = 0;
			double closestDistanceSquared = calculateDistanceSquared(means[0]);
			for (size_t index = 1; index < means.size(); index++) {
				// Assume all means have the same dimension as this point
				double distanceSquared = calculateDistanceSquared(means[index]);
				if (distanceSquared < closestDistanceSquared) {
					closestIndex = index;
					closestDistanceSquared = distanceSquared;
				}
			}

			return closestIndex;
		}

		// Returns 1 if A is closer, 2 is B is closer, and 0 if neither is.
		int chooseClosestPoint(const Point& a, const Point& b) {
			double distanceToA = calculateDistanceSquared(a);
			double distanceToB = calculateDistanceSquared(b);
			return (distanceToA < distanceToB) ? 1 : (distanceToB < distanceToA ? 2 : 0);
		}
};

class Bounds {
	private:
		double *min;
		double *max;
		int dimensionCount;

	public:
		Bounds(int dimensionCount): dimensionCount(dimensionCount) {}

		Bounds(const Bounds& bounds) {
			min = new double[bounds.dimensionCount];
			max = new double[bounds.dimensionCount];

			for (int i = 0; i < bounds.dimensionCount; i++) {
				min[i] = bounds.min[i];
				max[i] = bounds.max[i];
			}

			dimensionCount = bounds.dimensionCount;
		}

		double minimum(int dimension) const {
			if (dimension < 0 || dimension >= dimensionCount) {
				throw std::string("Invalid axis for bounds.min[]");
			} else {
				return min[dimension];
			}
		}

		double maximum(int dimension) const {
			if (dimension < 0 || dimension >= dimensionCount) {
				throw std::string("Invalid axis for bounds.max[]");
			} else {
				return max[dimension];
			}
		}

		// Using the perpendicular bisector of the line segment connecting A and B,
		// see if this rectangle is fully on one side of the bisector or the other side.
		int whichPointDominates(const Point& a, const Point& b) {
			// Each corner gets classified, and then we compare the
			// classifications of the corners of these boundaries.

			// A bitstring, one bit per dimension.
			// 0 = Minimum side, 1 = Maximum side.
			int cornerNumber = 0;
			double *cornerDimensions = new double[dimensionCount];
			Point corner = Point(cornerDimensions, dimensionCount);

			int dominantPoint = 0;

			// Iterate over the corners.
			for (cornerNumber = 0; cornerNumber < (1 << dimensionCount); cornerNumber++) {
				for (int i = 0; i < dimensionCount; i++) {
					bool isUpperBoundary = (cornerNumber & (1 << i));
					cornerDimensions[i] = isUpperBoundary ? max[i] : min[i];
				}
				int closestPoint = corner.chooseClosestPoint(a, b);
				if (dominantPoint == 0) {
					// Initialize the dominant point
					dominantPoint = closestPoint;
				} else if (dominantPoint != closestPoint) {
					// There is no dominant point
					return 0;
				}
			}

			return dominantPoint;
		}

		Bounds createClippedVersionAfterPoint(const Point& clipper, int axis) const {
			if (axis < 0 || axis >= dimensionCount) {
				throw std::string("Invalid axis for clipped version of bounds");
			} else {
				Bounds result = Bounds(*this);
				result.min[axis] = clipper[axis];
				return result;
			}
		}

		Bounds createClippedVersionBeforePoint(const Point& clipper, int axis) const {
			if (axis < 0 || axis >= dimensionCount) {
				throw std::string("Invalid axis for clipped version of bounds");
			} else {
				Bounds result = Bounds(*this);
				result.max[axis] = clipper[axis];
				return result;
			}
		}
		
		static Bounds defaultBounds(int dimensionCount) {
			auto result = Bounds(dimensionCount);
			result.min = new double[dimensionCount];
			result.max = new double[dimensionCount];
			for (int i = 0; i < dimensionCount; i++) {
				result.min[i] = 0;
				result.max[i] = 1;
			}
			return result;
		}
};

class PPM {
	private:
		Pixel **pixels;
		int width, height;

	public:
		PPM() {}

		PPM(int width, int height): width(width), height(height) {
			pixels = new Pixel*[height];
			for (int y = 0; y < height; y++) {
				pixels[y] = new Pixel[width];
				for (int x = 0; x < width; x++) {
					pixels[y][x] = white;
				}
			}
		}

		void setPixel(int x, int y, Pixel value) {
			if (x < 0 || x >= width || y < 0 || y >= height) {
				// std::cerr << "ERROR: Invalid Pixel " << x << ", " << y << "\n";
			} else {
				pixels[y][x] = value;
			}
		}

		Pixel getPixel(int x, int y) const {
			return pixels[y][x];
		}

		void drawCircle(int center_x, int center_y, double radius, Pixel color) {
			// starts with the topmost point
			int x = 0;
			int y = (int)(radius + 0.5); // round up

			int y2 = y * y;
			int y2_new = y2;

			int two_y = 2 * y - 1;

			while (y >= x) {
				// when X increases, see if the Y value should decrease.
				if ((y2 - y2_new) >= two_y) {
					y2 -= two_y;

					// decrease Y and 2Y
					y -= 1;
					two_y -= 2;
				}

				setPixel(x + center_x, y + center_y, color);
				setPixel(x + center_x, -y + center_y, color);
				setPixel(-x + center_x, y + center_y, color);
				setPixel(-x + center_x, -y + center_y, color);

				setPixel(y + center_x, x + center_y, color);
				setPixel(y + center_x, -x + center_y, color);
				setPixel(-y + center_x, x + center_y, color);
				setPixel(-y + center_x, -x + center_y, color);

				y2_new -= 2 * x - 3;
				x += 1;
			}
		}

		int getWidth() const { return width; }
		int getHeight() const { return height; }

		void getAllPixelColors(std::vector<Point>& colors) const {
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					const Pixel pixel = pixels[y][x];
					colors.push_back({(double)pixel.r, (double)pixel.g, (double)pixel.b});
				}
			}
		}

		void drawHorizontalLine(int y, int minimumX, int maximumX) {
			for (int x = minimumX; x < maximumX; x++) {
				setPixel(x, y, red);
			}
		}

		void drawVerticalLine(int x, int minimumY, int maximumY) {
			for (int y = minimumY; y < maximumY; y++) {
				setPixel(x, y, blue);
			}
		}

		friend std::ostream& operator<<(std::ostream& stream, const PPM& ppm);
		friend std::istream& operator>>(std::istream& stream, PPM& ppm);
};

std::istream& operator>>(std::istream& stream, PPM& ppm) {
	stream.ignore(2); // "P3"
	int maxActivation = 0, r, g, b;
	stream >> ppm.width >> ppm.height >> maxActivation;
	ppm.pixels = new Pixel*[ppm.height];
	for (int y = 0; y < ppm.height; y++) {
		ppm.pixels[y] = new Pixel[ppm.width];
		for (int x = 0; x < ppm.width; x++) {
			stream >> r >> g >> b;
			ppm.pixels[y][x] = {r, g, b};
		}
	}

	return stream;
}

std::ostream& operator<<(std::ostream& stream, const PPM& ppm) {
	// PPM header
	stream << "P3 ";
	// Image dimensions
	stream << ppm.width << ' ' << ppm.height << ' ';
	// Max value of each rgb channel
	stream << 255 << ' ';
	
	// Write the pixels like this:
	// R G B R G B R G B
	// R G B R G B R G B
	// Each R G B corresponds to a single pixel
	for (int y = 0; y < ppm.height; y++) {
		for (int x = 0; x < ppm.width; x++) {
			Pixel p = ppm.getPixel(x, y);
			stream << p.r << ' ' << p.g << ' ' << p.b << ' ';
		}
		stream << '\n';
	}

	return stream;
}

int cycleForward(int x, int max) {
	return (x + 1) % max;
}

class KDTreeNode {
	private:
		const int axis;
		const int dimensionCount;
		KDTreeNode *gt = nullptr, *lt = nullptr;
		Point value;
		Bounds boundaries;
		bool hasValue = false;

		KDTreeNode(int dimensionCount, int axis, Bounds boundaries, Point &value):
			dimensionCount(dimensionCount),
			axis(axis),
			value(value),
			boundaries(boundaries),
			hasValue(true) {}

	public:
		KDTreeNode(int dimensionCount):
			dimensionCount(dimensionCount),
			axis(0),
			boundaries(Bounds::defaultBounds(dimensionCount)) {}

		bool hasPoint() const {
			return hasValue;
		}

		const Point getPoint() const {
			return this->value;
		}

		int getAxis() const {
			return this->axis;
		}

		int size() const {
			int s = hasValue;
			if (gt != nullptr) s += gt->size();
			if (lt != nullptr) s += lt->size();
			return s;
		}

		void addPoint(Point &pointToAdd) {
			if (!hasValue) {
				this->value = pointToAdd;
				hasValue = true;
			} else {
				bool isGreaterThanBreakpoint = pointToAdd[axis] > this->value[axis];
				if (isGreaterThanBreakpoint) {
					if (gt == nullptr) {
						int nextAxis = cycleForward(this->axis, pointToAdd.getDimensionCount());

						gt = new KDTreeNode(
							dimensionCount,
							nextAxis,
							boundaries.createClippedVersionAfterPoint(pointToAdd, nextAxis),
							pointToAdd);

					} else {
						gt->addPoint(pointToAdd);
					}
				} else {
					if (lt == nullptr) {
						int nextAxis = cycleForward(this->axis, pointToAdd.getDimensionCount());
						
						lt = new KDTreeNode(
							dimensionCount,
							nextAxis,
							boundaries.createClippedVersionBeforePoint(pointToAdd, nextAxis),
							pointToAdd);

					} else {
						lt->addPoint(pointToAdd);
					}
				}
			}
		}

		void drawToPPM(PPM& ppm) {
			if (hasValue) {
				bool isHorizontal = axis == 0;
				if (isHorizontal) {
					ppm.drawHorizontalLine((int) (value[1] * 800), boundaries.minimum(0), boundaries.maximum(0));
				} else {
					ppm.drawVerticalLine((int) (value[0] * 800), boundaries.minimum(1), boundaries.maximum(1));
				}

				if (gt != nullptr) {
					gt->drawToPPM(ppm);
				}

				if (lt != nullptr) {
					lt->drawToPPM(ppm);
				}

				ppm.drawCircle((int)(value[0] * 800), (int)(value[1] * 800), 2, black);
			}
		}

		void printTree(int indent = 0) {
			std::string space(indent, ' ');
			std::cout << space << "K-D Tree along axis " << axis;
			if (hasValue) {
				std::cout << " with point (" << value[0] << ", " << value[1] << ")";
			}
			std::cout << '\n';

			if (gt != nullptr) {
				std::cout << space << "GTTree:\n";
				gt->printTree(indent + 2);
			}

			if (lt != nullptr) {
				std::cout << space << "LTTree:\n";
				lt->printTree(indent + 2);
			}
		}
};

std::ostream& operator<<(std::ostream& stream, const Point& point) {
	stream << '(';
	int dims = point.getDimensionCount();
	for (int i = 0; i < dims; i++) {
		stream << point[i];
		if (i < dims - 1) {
			stream << ", ";
		}
	}

	stream << ')';
	return stream;
}

double randomDouble() {
	return rand() / (double)(RAND_MAX);
}

void addRandomPoints(Points& out, int n) {
	for (int i = 0; i < n; i++) {
		out.push_back(Point(randomDouble(), randomDouble()));
	}
}

void classifyPoints(
	std::vector<Cluster>& clusters,
	const Points& points,
	const Means& means
) {
	for (Point point : points) {
		int classification = point.classify(means);
		clusters[classification].push_back(point);
	}
}

Point calcMean(const Cluster& cluster) {
	double sumX = 0;
	double sumY = 0;
	for (Point p : cluster) {
		sumX += p[0];
		sumY += p[1];
	}
	return Point(sumX / cluster.size(), sumY / cluster.size());
}

void calcMeans(Means& out, const std::vector<Cluster>& clusters) {
	for (Cluster cluster : clusters) {
		out.push_back(calcMean(cluster));
	}
}

#include <iomanip>
#include <fstream>
void savePoints(const std::vector<Point>& points, const char* filename = "points.txt") {
	std::ofstream out(filename);
	
	// makes all numbers write with exactly 17 digits after the decimal point
	out << std::fixed << std::setprecision(17);
	for (auto it = points.begin(); it != points.end(); it++) {
		out << it->at(0) << "  " << it->at(1) << "\n";
	}

	out.close();
}

void generateAndSavePoints(int number = 10) {
	std::vector<Point> points;
	addRandomPoints(points, number);
	savePoints(points);
}

void readPoints(std::vector<Point>& output, const char* filename = "points.txt") {
	std::ifstream in(filename);
	double x, y;
	while (in >> x) {
		in >> y;
		output.push_back(Point(x, y));
	}
}

void part4() {
	std::cout << "Lab 4 Part 4 - KD Trees and K Means!\n";
	std::cout << "Would you like to generate a new file of 10 points?\n";
	std::cout << "(Warning: Previous file will be overwritten)\n";
	std::cout << "Choose (yes/no): ";
	std::string answer;
	std::cin >> answer;
	if (answer == "yes") {
		generateAndSavePoints(10);
	}

	Points points;
	readPoints(points);

	KDTreeNode tree(2);

	for (Point& point : points) {
		tree.addPoint(point);
	}

	// tree.printTree();

	PPM ppm(800, 800);
	tree.drawToPPM(ppm);

	std::ofstream out("diagram.ppm");
	out << ppm;
	out.close();
}

#include <time.h>
int main() {
	std::srand(time(NULL));

	part4();	
}
