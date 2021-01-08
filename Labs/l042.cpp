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
			for (size_t i = 0; i < dimensions; i++) {
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
			for (size_t meanIndex = 1; meanIndex < means.size(); meanIndex++) {
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

typedef struct {
	int r, g, b;
} Pixel;

const Pixel black  = { 0x00, 0x00, 0x00 };
const Pixel white  = { 0xFF, 0xFF, 0xFF };
const Pixel red    = { 0xFF, 0x00, 0x00 };
const Pixel orange = { 0xFF, 0x7F, 0x00 };
const Pixel yellow = { 0xFF, 0x00, 0xFF };
const Pixel green  = { 0x00, 0xFF, 0x00 };
const Pixel blue   = { 0x00, 0x00, 0xFF };
const Pixel purple = { 0xFF, 0x00, 0xFF };

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

double getRandom() {
	return rand() / (double)(RAND_MAX);
}

void addRandomPoints(Points& out, int n) {
	for (int i = 0; i < n; i++) {
		out.push_back(Point(getRandom(), getRandom()));
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

#include <fstream>
void part2() {
	Points points;

	PPM toConvert;
	std::ifstream toConvertFile("./peppers.ppm");

	toConvertFile >> toConvert;
	toConvert.getAllPixelColors(points);

	std::cout << "Num points: " << points.size() << "\n";

	const int k = 4;
	Means means, meansTemporary;
	means.push_back(Point(  0,  0,  0));
	means.push_back(Point( 85, 85, 85));
	means.push_back(Point(170,170,170));
	means.push_back(Point(255,255,255));

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

	std::cout << "Found means\n";

	PPM ppm(toConvert.getWidth(), toConvert.getHeight());
	for (int x = 0; x < toConvert.getWidth(); x++) {
		for (int y = 0; y < toConvert.getHeight(); y++) {
			const Pixel pixel = toConvert.getPixel(x,y);

			const Point p = {
				(double) pixel.r,
				(double) pixel.g,
				(double) pixel.b
			};

			int meanIndex = p.classify(means);

			const Point classification = means[meanIndex];

			ppm.setPixel(x, y, {
				(int)classification[0],
				(int)classification[1],
				(int)classification[2]
			});
		}
	}

	std::cout << "Converted file\n";

	std::ofstream kmeansOutput("peppers.out.ppm");
	kmeansOutput << ppm;
	kmeansOutput.close();

	std::cout << "Done\n";
}

#include <time.h>
int main() {
	std::srand(time(NULL));

	part2();	
}
