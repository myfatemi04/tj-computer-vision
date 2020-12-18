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

std::ostream& operator<<(std::ostream& stream, const Point& point) {
	stream << '(';
	int dims = point.getDimensions();
	for (int i = 0; i < dims; i++) {
		stream << point[i];
		if (i < dims - 1) {
			stream << ", ";
		}
	}

	stream << ')';
	return stream;
}

typedef const int* pixel;

const pixel black  = new int[3] { 0x00, 0x00, 0x00 };
const pixel white  = new int[3] { 0xFF, 0xFF, 0xFF };
const pixel red    = new int[3] { 0xFF, 0x00, 0x00 };
const pixel orange = new int[3] { 0xFF, 0x7F, 0x00 };
const pixel yellow = new int[3] { 0xFF, 0x00, 0xFF };
const pixel green  = new int[3] { 0x00, 0xFF, 0x00 };
const pixel blue   = new int[3] { 0x00, 0x00, 0xFF };
const pixel purple = new int[3] { 0xFF, 0x00, 0xFF };

class PPM {
	private:
		pixel **pixels;
		int width, height;

	public:
		PPM(int width, int height): width(width), height(height) {
			pixels = new pixel*[height];
			for (int y = 0; y < height; y++) {
				pixels[y] = new pixel[width];
				for (int x = 0; x < width; x++) {
					pixels[y][x] = white;
				}
			}
		}

		void setPixel(int x, int y, pixel value) {
			if (x < 0 || x >= width || y < 0 || y >= height) {
				std::cerr << "ERROR: Invalid Pixel " << x << ", " << y << "\n";
			}
			pixels[y][x] = value;
		}

		pixel getPixel(int x, int y) const {
			return pixels[y][x];
		}

		int getWidth() const { return width; }
		int getHeight() const { return height; }

		friend std::ostream& operator<<(std::ostream& stream, const PPM& ppm);
};

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
			pixel p = ppm.getPixel(x, y);
			stream << p[0] << ' ' << p[1] << ' ' << p[2] << ' ';
		}
		stream << '\n';
	}

	return stream;
}

void bruteForce(const Points& points, int k) {
}

double random() {
	return rand() / (double)(RAND_MAX + 1);
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
#include <fstream>

int main() {
	std::srand(time(NULL));

	Points points;
	addRandomPoints(points, 50000);

	const int k = 5;
	Means means;
	Means meansTemporary;
	addRandomPoints(means, k);

	std::vector<Cluster> clusters;
	for (int i = 0; i < k; i++) {
		clusters.push_back(Cluster());
	}
	while (true) {
		for (int i = 0; i < k; i++) {
			clusters[i].clear();
		}
		meansTemporary.clear();

		classifyPoints(clusters, points, means);
		calcMeans(meansTemporary, clusters);

		if (means != meansTemporary) {
			means = meansTemporary;
		} else {
			break;
		}
	}

	for (Point p : means) {
		std::cout << p << '\n';	
	}

	PPM ppm(800, 800);

	const pixel *clusterColors = new pixel[5] {
		red,
		green,
		blue,
		orange,
		purple
	};

	for (int i = 0; i < k; i++) {
		Cluster cluster = clusters[i];
		std::cout << cluster.size() << ' ';
		for (Point point : cluster) {
			ppm.setPixel((int)(point[0] * 800), (int)(point[1] * 800), clusterColors[i]);
		}
	}

	std::ofstream kmeansOutput("kmeansOutput.ppm");
	kmeansOutput << ppm;
	kmeansOutput.close();

	std::cout << '\n';
	std::cout << "Done\n";
}
