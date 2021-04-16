// l053.cpp
// Michael Fatemi

#include <iostream>
#include <fstream>
#include <cmath>
#include <set>
#include <vector>

#define USE_DEBUG true
#if USE_DEBUG
# define dbg(x) std::cout << x
#else
# define dbg(x)
#endif

int max(int a, int b) {
	return (a > b) ? a : b;
}

int min(int a, int b) {
	return (a < b) ? a : b;
}
namespace tjcv {
	class GrayscaleImage;
	class ColorImage;

	class ColorImage {
		private:
			int*** pixels;
			int width, height;
			int max = 255;
		public:
			ColorImage(int width, int height, int max = 255);

			void set(int x, int y, int *value);
			int *get(int x, int y);

			int getWidth() const;
			int getHeight() const;
			std::pair<int, int> getSize() const;

			void save(std::string filename);

			GrayscaleImage toGrayscale();
			ColorImage clone();

			static ColorImage fromPPM(std::string filename);
	};

	class GrayscaleImage {
		private:
			int** pixels;
			int width, height;
			int max = 255;

		public:
			GrayscaleImage(int **pixels, int width, int height, int max = 255);
			GrayscaleImage(int width, int height, int max = 255);
	
			void save(std::string filename);
			
			void set(int x, int y, int value);
			int get(int x, int y);

			int getMax();
			int getWidth();
			int getHeight();
			std::pair<int, int> getSize();

			ColorImage toColor();

			GrayscaleImage convolve(GrayscaleImage filter);
	};

	ColorImage::ColorImage(int width, int height, int max): width(width), height(height), max(max) {
		this->pixels = new int**[height];
		for (int y = 0; y < height; y++) {
			this->pixels[y] = new int*[width];
			for (int x = 0; x < width; x++) {
				// Set to color white
				this->pixels[y][x] = new int[3] { 0, 0, 0 };
			}
		}
	}

	void ColorImage::set(int x, int y, int *value) {
		if (x < 0 || y < 0) return;
		if (x >= width || y >= height) return;
		this->pixels[y][x] = value;
	}

	int *ColorImage::get(int x, int y) {
		if (x < 0 || y < 0) return nullptr;
		if (x >= width || y >= height) return nullptr;
		return this->pixels[y][x];
	}

	int ColorImage::getWidth() const { return this->width; }
	int ColorImage::getHeight() const { return this->height; }
	std::pair<int, int> ColorImage::getSize() const { return { this->width, this->height }; }

	void ColorImage::save(std::string filename) {
		std::ofstream handle(filename);
		// P3 [width] [height] [max intensity]
		handle << "P3 " << width << " " << height << " 255\n";
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				auto pixel = pixels[y][x];
				int r = pixel[0];
				int g = pixel[1];
				int b = pixel[2];
				handle << r << " " << g << " " << b << " ";
			}
			handle << '\n';
		}

		handle.close();
	}

	GrayscaleImage ColorImage::toGrayscale() {
		GrayscaleImage result(width, height, max);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				result.set(x, y, (pixels[y][x][0] + pixels[y][x][1] + pixels[y][x][2]) / 3);
			}
		}
		return result;
	}

	ColorImage ColorImage::clone() {
		ColorImage cloned(width, height);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				auto pixel = get(x, y);
				cloned.set(x, y, new int[3] {
					pixel[0],
					pixel[1],
					pixel[2]
				});
			}
		}
		
		return cloned;
	}

	ColorImage ColorImage::fromPPM(std::string filename) {
		std::ifstream handle(filename);
		std::string _ppmtype;
		handle >> _ppmtype;
		int width, height, _max;
		handle >> width >> height >> _max;

		ColorImage image(width, height, _max);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int *pixel = new int[3];
				handle >> pixel[0] >> pixel[1] >> pixel[2];
				image.set(x, y, pixel);
			}
		}
		handle.close();
		return image;
	}

	GrayscaleImage::GrayscaleImage(int **pixels, int width, int height, int max):
		pixels(pixels),
		width(width),
		height(height),
		max(max) {}

	GrayscaleImage::GrayscaleImage(int width, int height, int max): width(width), height(height), max(max) {
		this->pixels = new int*[height];
		for (int y = 0; y < height; y++) {
			this->pixels[y] = new int[width];
			for (int x = 0; x < width; x++) {
				this->pixels[y][x] = 0;
			}
		}
	}
	
	void GrayscaleImage::save(std::string filename) {
		std::ofstream handle(filename);
		// P2 [width] [height] [max intensity]
		handle << "P2 " << width << " " << height << ' ' << max << '\n';
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				handle << abs(pixels[y][x]) << " ";
			}
			handle << '\n';
		}

		handle.close();
	}
	
	void GrayscaleImage::set(int x, int y, int value) {
		if (x < 0 || y < 0) return;
		if (x >= width || y >= height) return;
		this->pixels[y][x] = value;
	}

	int GrayscaleImage::get(int x, int y) {
		if (x < 0 || y < 0) return -1;
		if (x >= width || y >= height) return -1;
		return this->pixels[y][x];
	}

	int GrayscaleImage::getMax() { return this->max; }
	int GrayscaleImage::getWidth() { return this->width; }
	int GrayscaleImage::getHeight() { return this->height; }
	std::pair<int, int> GrayscaleImage::getSize() { return { this->width, this->height }; }

	ColorImage GrayscaleImage::toColor() {
		ColorImage color(width, height);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int intensity = pixels[y][x];
				// R, G, and B are all the same value
				color.set(x, y, new int[3] { intensity, intensity, intensity });
			}
		}

		return color;
	}

	GrayscaleImage GrayscaleImage::convolve(GrayscaleImage filter) {
		dbg("Convolving an image around a filter " << filter.getWidth() << "x" << filter.getHeight() << '\n');

		GrayscaleImage convolved(width, height, max);
		int filterSum = 0;
		for (int i = 0; i < filter.getHeight(); i++) {
			for (int j = 0; j < filter.getWidth(); j++) {
				filterSum += filter.get(i, j);
			}
		}
		dbg("Filter sum: " << filterSum << '\n');

		int filterRadiusX = filter.getWidth() >> 1;
		int filterRadiusY = filter.getHeight() >> 1;
		dbg("Filter radius: " << filterRadiusX << ", " << filterRadiusY << '\n');

		// Initialize the edges to 0
		for (int y = 0; y < height; y++) {
			for (int i = 0; i < filterRadiusX; i++) {
				convolved.set(i, y, 0);
				convolved.set(width - i - 1, y, 0);
			}
		}

		for (int x = 0; x < width; x++) {
			for (int i = 0; i < filterRadiusY; i++) {
				convolved.set(x, i, 0);
				convolved.set(x, height - i - 1, 0);
			}
		}

		for (int y = filterRadiusY; y < height - filterRadiusY; y++) {
			for (int x = filterRadiusX; x < width - filterRadiusX; x++) {
				int total = 0;
				for (int relativeX = -filterRadiusX; relativeX <= filterRadiusX; relativeX++) {
					for (int relativeY = -filterRadiusY; relativeY <= filterRadiusY; relativeY++) {
						total += get(y + relativeY, x + relativeX) * filter.get(relativeY + filterRadiusY, relativeX + filterRadiusX);
					}
				}

				convolved.set(x, y, filterSum != 0 ? (total / filterSum) : total);
			}
		}
		
		return convolved;
	}
	typedef struct {
		int x, y;
		double radius;
	} Circle;

	class BresenhamPixelIterator {
		private:
			int startX, startY;
			int currentX, currentY;
			double buildup;
			double angle;
			bool iterateOverX;
			
		public:
			BresenhamPixelIterator(int startX, int startY, double angle): startX(startX), startY(startY), angle(angle) {
				currentX = startX;
				currentY = startY;
				buildup = 0;
				
				iterateOverX = abs(cos(angle)) > abs(sin(angle));
			}

			void reset() {
				currentX = startX;
				currentY = startY;
				buildup = 0;
			}

			int getX() const {
				return currentX;
			}

			int getY() const {
				return currentY;
			}
			
			void step(int count) {
				if (iterateOverX) {
					currentX += count;
					buildup += count * sin(angle) / cos(angle);
					int overshoot = (int) buildup;
					if (overshoot != 0) {
						buildup -= overshoot;
						currentY += overshoot;
					}
				} else {
					currentY += count;
					buildup += count * sin(angle) / cos(angle);
					int overshoot = (int) buildup;
					if (overshoot != 0) {
						buildup -= overshoot;
						currentX += overshoot;
					}
				}
			}
		};

	/**
	 * getCirclePixels returns a vector of (x, y) pixels. It uses an iterative, symmetric method of
	 * finding the pixels of a circle. The center must be provided as a pair of integers, but the
	 * radius may be a decimal value.
	 */
	std::vector<std::pair<int, int>> getCirclePixels(int centerX, int centerY, double radius) {
		// starts with the topmost point
		int x = 0;
		int y = round(radius);

		int y2 = y * y;
		int y2_new = y2;

		int two_y = 2 * y - 1;

		std::vector<std::pair<int, int>> pixels;

		while (y >= x) {
			// when X increases, see if the Y value should decrease.
			if ((y2 - y2_new) >= two_y) {
				y2 -= two_y;

				// decrease Y and 2Y
				y -= 1;
				two_y -= 2;
			}

			pixels.push_back({  x + centerX,  y + centerY });
			pixels.push_back({  x + centerX, -y + centerY });
			pixels.push_back({ -x + centerX,  y + centerY });
			pixels.push_back({ -x + centerX, -y + centerY });

			pixels.push_back({  y + centerX,  x + centerY });
			pixels.push_back({  y + centerX, -x + centerY });
			pixels.push_back({ -y + centerX,  x + centerY });
			pixels.push_back({ -y + centerX, -x + centerY });

			y2_new -= 2 * x - 3;

			x += 1;
		}

		return pixels;
	}

	double **make2DDoubleArray(int dim1, int dim2) {
		double **arr = new double*[dim1];
		for (int i = 0; i < dim1; i++) {
			arr[i] = new double[dim2];
		}

		return arr;
	}
	
	int **make2DIntArray(int dim1, int dim2) {
		int **arr = new int*[dim1];
		for (int i = 0; i < dim1; i++) {
			arr[i] = new int[dim2];
			for (int j = 0; j < dim2; j++) {
				arr[i][j] = 0;
			}
		}

		return arr;
	}

	bool inbounds(int x, int y, int maxX, int maxY) {
		return (x >= 0 && x < maxX) && (y >= 0 && y < maxY);
	}

	const double SIN45 = sqrt(2) / 2;

	int countPixelsToDrawCircle(int radius) {
		// Our circle rasterization method uses symmetry
		// There is one pixel on every X value for the first 45 degrees of the circle
		// (from 12:00 to 1:30 on a clock.)
		// So, we take the X length of the first 45 degrees of the circle, and multiply it by 8
		// to obtain the number of pixels needed to draw the circle.

		return (int) (SIN45 * radius * 8);
	}

	void drawCircle(ColorImage image, int x, int y, double radius, int* color) {
		auto pixels = getCirclePixels(x, y, radius);
		for (const auto& pixel : pixels) {
			int px = pixel.first;
			int py = pixel.second;
			image.set(px, py, color);
		}
	}

	void drawFilledCircle(ColorImage image, int x, int y, int radius, int* color) {
		for (int r = 0; r < radius; r++) {
			drawCircle(image, x, y, r, color);
		}
	}
}

namespace lab5 {
	using tjcv::GrayscaleImage;
	using tjcv::ColorImage;
	typedef struct {
		GrayscaleImage edges;
		double** angles;
		GrayscaleImage xGradient, yGradient;
	} EdgeDetectionResult;

	GrayscaleImage verticalSobel {
		new int*[3] {
			new int[3] {	1,	0, -1 },
			new int[3] {	2,	0, -2 },
			new int[3] {	1,	0, -1 }
		},
		3,
		3
	}, horizontalSobel {
		new int*[3] {
			new int[3] {	1,	2,	1 },
			new int[3] {	0,	0,	0 },
			new int[3] { -1, -2, -1 }
		},
		3,
		3
	}, gaussian3 {
		new int*[3] {
			new int[3] { 1, 2, 1 },
			new int[3] { 2, 4, 2 },
			new int[3] { 1, 2, 1 },
		},
		3,
		3
	}, gaussian5 {
		new int*[5] {
			new int[5] { 1,  4,  7,  4, 1 },
			new int[5] { 4, 16, 26, 16, 4 },
			new int[5] { 7, 26, 41, 26, 7 },
			new int[5] { 4, 16, 26, 16, 4 },
			new int[5] { 1,  4,  7,  4, 1 },
		},
		5,
		5
	};
	
	/**
	 * hysteresis() takes an image and runs a double threshold on it. Then, any "weak" edges
	 * touching a "strong" edge are automatically promoted to "strong" edges. This process is
	 * repeated until there are no more changes.
	 */
	GrayscaleImage hysteresis(GrayscaleImage magnitudes, int lowerThreshold, int upperThreshold) {
		// set of (x, y) pairs
		std::set<std::pair<int, int>> unvisited;
		GrayscaleImage newImage(magnitudes.getWidth(), magnitudes.getHeight(), 1);

		for (int y = 0; y < magnitudes.getHeight(); y++) {
			for (int x = 0; x < magnitudes.getWidth(); x++) {
				int pixelValue = magnitudes.get(x, y);
				if (pixelValue > upperThreshold) {
					unvisited.insert({ x, y });
					newImage.set(x, y, 2);
				} else if (pixelValue > lowerThreshold) {
					newImage.set(x, y, 1);
				} else {
					newImage.set(x, y, 0);
				}
			}
		}

		// return newImage;

		while (unvisited.size()) {
			// Pop a location from the set
			auto next = *unvisited.begin();
			int x = next.first, y = next.second;
			unvisited.erase(unvisited.begin());

			// Check nearby locations. x_ and y_ represent points around the current pixel.
			for (int x_ = x - 1; x_ <= x + 1; x_++) {
				// Boundary check
				if (x_ < 0 || x_ >= magnitudes.getWidth()) continue;
				for (int y_ = y - 1; y_ <= y + 1; y_++) {
					// Boundary check
					if (y_ < 0 || y_ >= magnitudes.getHeight()) continue;

					// If the current value is 1 (reached lower threshold, but not upper threshold),
					// then because it touches a strong edge with a value of 2, it is automatically
					// promoted. Then, we add it to the queue to update its neighbors as well.
					if (newImage.get(x_, y_) == 1) {
						newImage.set(x_, y_, 2);
						unvisited.insert({x_, y_});
					}
				}
			}
		}

		for (int y = 0; y < newImage.getHeight(); y++) {
			for (int x = 0; x < newImage.getWidth(); x++) {
				// Convert 2 to 1, and 1 to 0.
				newImage.set(x, y, newImage.get(x, y) >> 1);
			}
		}

		return newImage;
	}

	double findMagnitudeSquared(double a, double b) {
		return a * a + b * b;
	}

	bool inBounds(int x, int y, int width, int height) {
		return x > 0 && y > 0 && x < width && y < height;
	}

	/**
	 * nonMaxSuppression() looks at the magnitude and directions of each edge gradient.
	 * If a pixel along an edge is not the strongest edge along its gradient, it is suppressed.
	 */
	GrayscaleImage nonMaxSuppression(GrayscaleImage xGradient, GrayscaleImage yGradient, GrayscaleImage magnitudes) {
		int width = xGradient.getWidth();
		int height = xGradient.getHeight();
		GrayscaleImage im(width, height, 1);
		for (int y = 0; y < xGradient.getWidth(); y++) {
			for (int x = 0; x < xGradient.getHeight(); x++) {
				im.set(x, y, 0);
				double currentMagnitude = magnitudes.get(x, y);

				// angle is in the range [-pi, pi].
				double angleRadians = atan2(yGradient.get(x, y), xGradient.get(x, y));
				double angleDegrees = angleRadians / (3.1415926535897932) / 2 * 360;

				// Find the amount to go in the X and Y directions to follow the gradient

				int dx = 0, dy = 0;
				if (angleDegrees > 22.5 && angleDegrees < (180 - 22.5)) {
					dy = 1;
				} else if (angleDegrees < -22.5 && angleDegrees > (-180 + 22.5)) {
					dy = -1;
				}

				{
					double absDegrees = abs(angleDegrees);
					if (absDegrees > (90 + 22.5)) {
						dx = -1;
					} else if (absDegrees < (90 - 22.5)) {
						dx = 1;
					}
				}
				
				// std::cout << angleDegrees << ": " << dx << ", " << dy << '\n';

				if (inBounds(x + dx, y + dy, width, height)) {
					if (currentMagnitude > magnitudes.get(x + dx, y + dy)) {
						if (inBounds(x - dx, y - dy, width, height)) {
							if (currentMagnitude > magnitudes.get(x - dx, y - dy)) {
								im.set(x, y, 1);
							}
						}
					}
				}
			}
		}
		
		return im;
	}

	GrayscaleImage combineSobel(GrayscaleImage first, GrayscaleImage second) {
		GrayscaleImage combined(first.getWidth(), first.getHeight(), first.getMax());
		for (int y = 0; y < first.getHeight(); y++) {
			for (int x = 0; x < first.getWidth(); x++) {
				combined.set(x, y, (int) sqrt(findMagnitudeSquared(first.get(x, y), second.get(x, y))));
			}
		}
		
		return combined;
	}

	GrayscaleImage applyThreshold(GrayscaleImage image, double threshold) {
		GrayscaleImage result(image.getWidth(), image.getHeight(), image.getMax());
		for (int y = 0; y < image.getHeight(); y++) {
			for (int x = 0; x < image.getWidth(); x++) {
				result.set(x, y, image.get(x, y) >= threshold ? image.getMax() : 0);
			}
		}

		return result;
	}

	GrayscaleImage combineImages(GrayscaleImage afterHysteresis, GrayscaleImage afterNonMaxSuppression) {
		GrayscaleImage newImage(afterHysteresis.getWidth(), afterHysteresis.getHeight(), afterHysteresis.getMax());

		for (int y = 0; y < afterHysteresis.getHeight(); y++) {
			for (int x = 0; x < afterHysteresis.getWidth(); x++) {
				bool bothValid = afterHysteresis.get(x, y) != 0 && afterNonMaxSuppression.get(x, y) != 0;
				newImage.set(x, y, bothValid);
			}
		}

		return newImage;
	}

	/**
	 * Returns a 2D array of radian angle measures with the same dimensions as the input images.
	 */
	double **getEdgeAnglesFromGradients(GrayscaleImage xGradient, GrayscaleImage yGradient) {
		double **measures = new double*[xGradient.getHeight()];
		for (int y = 0; y < xGradient.getHeight(); y++) {
			measures[y] = new double[xGradient.getWidth()];
			for (int x = 0; x < xGradient.getWidth(); x++) {
				measures[y][x] = atan2(yGradient.get(x, y), xGradient.get(x, y));
			}
		}

		return measures;
	}

	/**
	 * detectEdges() detects the edges of an image.
	 * @param grayscale The grayscale input image
	 * @param lowerThreshold The lower threshold to use for hysteresis (Weak edge)
	 * @param upperThreshold The upper threshold to use for hysteresis (Strong edge)
	 */
	EdgeDetectionResult detectEdges(GrayscaleImage grayscale, int lowerThreshold, int upperThreshold) {
		// GrayscaleImage afterGaussian = convolve(grayscale, gaussian5);
		GrayscaleImage xGradient = grayscale.convolve(horizontalSobel);
		GrayscaleImage yGradient = grayscale.convolve(verticalSobel);
		
		GrayscaleImage magnitudes = combineSobel(xGradient, yGradient);

		GrayscaleImage afterHysteresis = hysteresis(magnitudes, lowerThreshold, upperThreshold);
		GrayscaleImage afterNonMaxSuppression = nonMaxSuppression(xGradient, yGradient, magnitudes);

		return EdgeDetectionResult {
			combineImages(afterHysteresis, afterNonMaxSuppression),
			getEdgeAnglesFromGradients(xGradient, yGradient),
			xGradient,
			yGradient
		};
	}
}

namespace lab6 {
	using tjcv::GrayscaleImage;

	/**
	 * Using Bresenham's algorithm, casts votes along a line. The rayWidth option allows you to cast lines that are wider than one pixel.
	 */
	void castVotesForOnePixel(
		int **votes,
		int x,
		int y,
		int width,
		int height,
		double angle,
		int rayWidth = 0) {
		using tjcv::inbounds;

		for (int offset = -rayWidth; offset <= rayWidth; offset++) {
			// Offset perpendicular to the line
			// cos -> sin
			// sin -> -cos
			int offsetX = (int) (sin(angle) * offset);
			int offsetY = (int) (-cos(angle) * offset);

			tjcv::BresenhamPixelIterator it(x, y, angle);
			int cx, cy;
			while (cx = it.getX(), cy = it.getY(), inbounds(cx, cy, width, height)) {
				votes[cy][cx]++;
				it.step(1);
			}

			it.reset();

			while (cx = it.getX(), cy = it.getY(), inbounds(cx, cy, width, height)) {
				votes[cy][cx]++;
				it.step(-1);
			}
		}

	}

	int **castVotes(GrayscaleImage edges, double **angles, int rayWidth = 0) {
		int **votes = tjcv::make2DIntArray(edges.getHeight(), edges.getWidth());

		for (int y = 0; y < edges.getHeight(); y++) {
			for (int x = 0; x < edges.getWidth(); x++) {
				if (edges.get(x, y)) {
					castVotesForOnePixel(votes, x, y, edges.getWidth(), edges.getHeight(), angles[y][x], rayWidth);
				}
			}
		}

		return votes;
	}

	/**
	 * findCenters() returns a vector of (x, y) points. This method takes in an int**, which counts up the votes,
	 * and a width and height, to specify the area within the int** to search. The optional threshold parameter
	 * can be used to specify the minimum number of votes required for a point to be classified as a center.
	 */
	std::vector<std::pair<int, int>> findCenters(int **votes, int width, int height, int threshold = 10) {
		std::vector<std::pair<int, int>> centers;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (votes[y][x] >= threshold) {
					centers.push_back({x, y});
				}
			}
		}
		return centers;
	}

	/**
	 * scoreCircleCandidate scores the likelihood of a radius being a circle.
	 * The way it does this is by iterating over the pixels that would compose a circle.
	 * If a pixel is an edge, and also connected to other edges, then it will receive a
	 * higher score than other edges.
	 */
	double scoreCircleCandidate(
		GrayscaleImage edges,
		int x,
		int y,
		int radius,
		double unconnectedScore,
		double connectedScore) {
		double score = 0;
		auto pixels = tjcv::getCirclePixels(x, y, radius);
		for (const auto& pixel : pixels) {
			int px = pixel.first;
			int py = pixel.second;
			if (tjcv::inbounds(px, py, edges.getWidth(), edges.getHeight())) {
				if (edges.get(px, py) > 0) {
					// Check if it has a neighboring edge
					bool hasNeighbor = false;
					for (int i = -1; i <= 1; i++) {
						for (int j = -1; j <= 1; j++) {
							if (edges.get(i, j) > 0) {
								hasNeighbor = true;
								break;
							}
						}
					}
					
					score += hasNeighbor * connectedScore + !hasNeighbor * unconnectedScore;
				}
			}
		}

		return score;
	}

	/**
	 * findRadii returns the radii that contain a certain number of edge pixels along their circumference. The radius is checked on the interval [minRadius, maxRadius). minRatio specifies the minimum ratio of empty edges to fillled edges.
	 */
	std::vector<int> findRadii(
		GrayscaleImage edges,
		int x,
		int y,
		int minRadius,
		int maxRadius,
		int ringWidth,
		double minRatio) {
		// Count the edges on each radius
		double *edgesByRadius = new double[maxRadius - minRadius];
		for (int radius = minRadius; radius < maxRadius; radius++) {
			edgesByRadius[radius - minRadius] = scoreCircleCandidate(edges, x, y, radius, 0.5, 1);
		}
		
		std::vector<int> radii;

		// Iterate over each possible radius, using the given ringWidth.
		for (int radius = minRadius; radius < maxRadius; radius++) {
			// Keep track of how many edges are found along the circle
			int ringEdgeCount    = edgesByRadius[radius - minRadius];
			// Keep track of how many edges could possibly be found along the circle
			int ringMaxEdgeCount = tjcv::countPixelsToDrawCircle(radius);

			// Ring radii must be in the interval requested
			int ringStartRadius  = max(minRadius, radius - ringWidth);
			int ringEndRadius    = min(radius + ringWidth, maxRadius - 1);

			for (int ringRadius = ringStartRadius; ringRadius <= ringEndRadius; ringRadius++) {
				// ringMaxEdgeCount += tjcv::countPixelsToDrawCircle(ringRadius);
				ringEdgeCount    += edgesByRadius[ringRadius - minRadius];
			}

			double ratio = ((double) ringEdgeCount / ringMaxEdgeCount);
			if (radius < 30 && radius > 10) {
				dbg("found/max (=): " << ringEdgeCount << "/" << ringMaxEdgeCount << "(" << ratio << ")\n");
			}

			if (ratio > minRatio) {
				radii.push_back(radius);
			}
		}
		return radii;
	}

	GrayscaleImage createVotesGraph(int **votes, int width, int height) {
		int maxIntensity = 1;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (votes[y][x] > maxIntensity) {
					maxIntensity = votes[y][x];
				}
			}
		}

		return GrayscaleImage { votes, width, height, maxIntensity };
	}

	void part1() {
		using tjcv::ColorImage;
		dbg("Loading image\n");
		auto color = ColorImage::fromPPM("image.ppm");
		auto grayscale = color.toGrayscale();

		dbg("Detecting edges\n");
		auto detection = lab5::detectEdges(grayscale, 160, 180);
		detection.edges.save("imagef.ppm");

		dbg("Casting votes\n");
		auto votes = lab6::castVotes(detection.edges, detection.angles, 2);
		auto votesGraph = lab6::createVotesGraph(votes, detection.edges.getWidth(), detection.edges.getHeight());
		votesGraph.save("imagev.ppm");

		dbg("Finding centers\n");
		auto centers = lab6::findCenters(votes, grayscale.getWidth(), grayscale.getHeight(), 200);

		dbg("Found " << centers.size() << " centers\n");

		dbg("Drawing centers\n");
		auto colorWithCenters = color.clone();
		int* CENTER_COLOR = new int[3] { 255, 0, 0 };
		// On all the centers found, create a filled red circle of radius 5.
		for (const auto& center : centers) {
			int x = center.first;
			int y = center.second;
			// dbg("Center location: " << x << ", " << y << '\n');
			// Iterate over radii
			for (int r = 1; r < 5; r++) {
				for (const auto& pixel : tjcv::getCirclePixels(x, y, r)) {
					int px = pixel.first;
					int py = pixel.second;
					colorWithCenters.set(px, py, CENTER_COLOR);
				}
			}
		}

		dbg("Saving centers\n");
		colorWithCenters.save("imageCC.ppm");
	}

	void part2() {
		using tjcv::ColorImage;

		dbg("Loading image\n");
		auto colorImage = ColorImage::fromPPM("image.ppm");
		auto grayscaleImage = colorImage.toGrayscale();

		dbg("Detecting edges\n");
		auto detection = lab5::detectEdges(grayscaleImage, 160, 180);
		detection.edges.save("imagef.ppm");

		dbg("Casting votes\n");
		auto votes = lab6::castVotes(detection.edges, detection.angles, 2);
		auto votesGraph = lab6::createVotesGraph(votes, detection.edges.getWidth(), detection.edges.getHeight());
		votesGraph.save("imagev.ppm");

		dbg("Finding centers\n");
		auto centers = lab6::findCenters(votes, grayscaleImage.getWidth(), grayscaleImage.getHeight(), 200);

		dbg("Found " << centers.size() << " centers\n");

		int* CIRCLE_COLOR = new int[3] { 0, 255, 0 };
		int* CENTER_COLOR = new int[3] { 255, 0, 0 };

		dbg("Finding radii\n");
		int foundRadiusCount = 0;
		for (int i = 0; i < centers.size(); i++) {
			const auto& center = centers.at(i);
			int x = center.first;
			int y = center.second;
			// tjcv::drawFilledCircle(colorImage, x, y, 5, CENTER_COLOR);
			
			auto radii = lab6::findRadii(detection.edges, x, y, 10, 30, 2, 0.9);
			for (int radius : radii) {
				dbg("Circle: {x=" << x << ", y=" << y << ", r=" << radius << "}\n");
				tjcv::drawCircle(colorImage, x, y, radius, CIRCLE_COLOR);
			}

			foundRadiusCount += radii.size();
		}

		dbg("Found " << foundRadiusCount << " radii\n");
		dbg("Saving\n");

		colorImage.save("imagecircles.ppm");
	}
}

int main() {
	lab6::part2();
}
