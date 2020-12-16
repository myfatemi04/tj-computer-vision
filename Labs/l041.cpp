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

class Point {
	private:
		size_t dimensions;
		double *vars;

	public:
		Point(double x, double y) {
			this -> dimensions = 2;
			this -> vars = new double[2] {x, y};
		}

		Point(double x, double y, double z) {
			this -> dimensions = 3;
			this -> vars = new double[3] {x, y, z};
		}

		double *getVars() const {
			return vars;
		}

		int getDimensions() const {
			return dimensions;
		}

		double distance2(const Point& other) const {
			double sum = 0, difference;
			for (size_t dimension = 0; dimension < max(dimensions, other.dimensions); dimension++) {
				difference = other[dimension] - at(dimension);
				sum += difference * difference;
			}
			return sum;
		}

		double at(size_t dimension) const {
			if (dimension < dimensions) {
				return vars[dimension];
			} else {
				return 0;
			}
		}

		double operator[](size_t dimension) const {
			return at(dimension);
		}

		bool operator==(const Point& other) const {
			if (dimensions != other.dimensions) {
				return false;
			}
			for (int i = 0; i < dimensions; i++) {
				if (at(i) != other[i]) {
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
			double closestDistance2 = distance2(means[0]);
			for (int meanIndex = 1; meanIndex < means.size(); meanIndex++) {
				// Assume all means have the same dimension as this point
				double distance2_ = distance2(means[meanIndex]);
				if (distance2_ < closestDistance2) {
					closestIndex = meanIndex;
					closestDistance2 = distance2_;
				}
			}

			return closestIndex;
		}
};

typedef unsigned char byte;

class PPM {
	private:
		byte **pixels;
		int width, height;

	public:
		PPM(int width, int height): width(width), height(height) {
			pixels = new byte*[height];
			for (int y = 0; y < height; y++) {
				pixels[y] = new byte[width];
				for (int x = 0; x < width; x++) {
					pixels[y][x] = 0;
				}
			}
		}

		void setPixel(int x, int y, byte value) const {
			pixels[y][x] = value;
		}

		byte getPixel(int x, int y) const {
			return pixels[y][x];
		}
};

void bruteForce(const Points& points, int k) {
}

double random() {
	return rand() / (double)RAND_MAX;
}

void addRandomPoints(Points& out, int n) {
	for (int i = 0; i < n; i++) {
		out.push_back(Point(random(), random()));
	}
}

/**
 * 
 */
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

#include <time.h>

int main() {
	std::srand(time(NULL));

	Points points;
	addRandomPoints(points, 500);

	Means means; 
	addRandomPoints(means, 3);

	std::vector<Cluster> clusters;
	classifyPoints(clusters, points, means);

	Means newMeans;
	calcMeans(newMeans, clusters);

	if (newMeans != means) {
		std::cout << "new means are different\n";
	} else {
		std::cout << "means are same\n";
	}
}
