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

namespace tjcv {
	typedef struct {
		int*** pixels;
		int width, height;
	} ColorImage;

	typedef struct {
		int** pixels;
		int width, height;
		int max = 255;
	} GrayscaleImage;

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

	void saveGrayscalePPM(std::string filename, GrayscaleImage image) {
		std::ofstream handle(filename);
		// P2 [width] [height] [max intensity]
		handle << "P2 " << image.width << " " << image.height << ' ' << image.max << '\n';
		for (int y = 0; y < image.height; y++) {
			for (int x = 0; x < image.width; x++) {
				handle << abs(image.pixels[y][x]) << " ";
			}
			handle << '\n';
		}

		handle.close();
	}

	ColorImage loadColorPPM(std::string filename) {
		std::ifstream handle(filename);
		std::string _ppmtype;
		handle >> _ppmtype;
		int width, height, _max;
		handle >> width >> height >> _max;
		auto pixels = new int**[height];
		for (int y = 0; y < height; y++) {
			pixels[y] = new int*[width];
			for (int x = 0; x < width; x++) {
				pixels[y][x] = new int[3];
				handle >> pixels[y][x][0] >> pixels[y][x][1] >> pixels[y][x][2];
			}
		}
		handle.close();
		return { pixels, width, height };
	}

	void saveColorPPM(std::string filename, ColorImage image) {
		std::ofstream handle(filename);
		// P3 [width] [height] [max intensity]
		handle << "P3 " << image.width << " " << image.height << " 255\n";
		for (int y = 0; y < image.height; y++) {
			for (int x = 0; x < image.width; x++) {
				auto pixel = image.pixels[y][x];
				int r = pixel[0];
				int g = pixel[1];
				int b = pixel[2];
				handle << r << " " << g << " " << b << " ";
			}
			handle << '\n';
		}

		handle.close();
	}

	GrayscaleImage convertToGrayscale(ColorImage image) {
		int** pixels = new int*[image.height];
		for (int y = 0; y < image.height; y++) {
			pixels[y] = new int[image.width];
			for (int x = 0; x < image.width; x++) {
				pixels[y][x] = (image.pixels[y][x][0] + image.pixels[y][x][1] + image.pixels[y][x][2]) / 3;
			}
		}
		return { pixels, image.width, image.height };
	}

	/**
	 * This simply creates a color image by taking the gray value and using it as the R, G, and B values.
	 */
	ColorImage convertToColor(GrayscaleImage gray) {
		ColorImage color;
		color.width = gray.width;
		color.height = gray.height;
		color.pixels = new int**[gray.height];
		for (int y = 0; y < gray.height; y++) {
			color.pixels[y] = new int*[gray.width];
			for (int x = 0; x < gray.width; x++) {
				int intensity = gray.pixels[y][x];
				// R, G, and B are all the same value
				color.pixels[y][x] = new int[3] {
					intensity,
					intensity,
					intensity
				};
			}
		}

		return color;
	}

	/**
	 * This method creates a new color image
	 */
	ColorImage cloneColorImage(ColorImage image) {
		ColorImage cloned;
		cloned.width = image.width;
		cloned.height = image.height;
		cloned.pixels = new int**[image.height];
		for (int y = 0; y < image.height; y++) {
			cloned.pixels[y] = new int*[image.width];
			for (int x = 0; x < image.width; x++) {
				auto pixel = image.pixels[y][x];
				cloned.pixels[y][x] = new int[3] {
					pixel[0],
					pixel[1],
					pixel[2]
				};
			}
		}
		
		return cloned;
	}

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

	GrayscaleImage convolve(GrayscaleImage image, GrayscaleImage filter) {
		dbg("Convolving an image around a filter " << filter.width << "x" << filter.height << '\n');

		auto pixels = new int*[image.height];
		
		int filterSum = 0;
		for (int i = 0; i < filter.height; i++) {
			for (int j = 0; j < filter.width; j++) {
				filterSum += filter.pixels[i][j];
			}
		}
		dbg("Filter sum: " << filterSum << '\n');

		int filterRadiusX = filter.width >> 1;
		int filterRadiusY = filter.height >> 1;
		dbg("Filter radius: " << filterRadiusX << ", " << filterRadiusY << '\n');

		// Initialize the edges to 0
		for (int y = 0; y < image.height; y++) {
			pixels[y] = new int[image.width];
			for (int i = 0; i < filterRadiusX; i++) {
				pixels[y][i] = 0;
				pixels[y][image.width - i - 1] = 0;
			}
		}
		for (int x = 0; x < image.width; x++) {
			for (int i = 0; i < filterRadiusY; i++) {
				pixels[i][x] = 0;
				pixels[image.height - i - 1][x] = 0;
			}
		}

		for (int y = filterRadiusY; y < image.height - filterRadiusY; y++) {
			for (int x = filterRadiusX; x < image.width - filterRadiusX; x++) {
				int total = 0;
				for (int relativeX = -filterRadiusX; relativeX <= filterRadiusX; relativeX++) {
					for (int relativeY = -filterRadiusY; relativeY <= filterRadiusY; relativeY++) {
						total += image.pixels[y + relativeY][x + relativeX] * filter.pixels[relativeY + filterRadiusY][relativeX + filterRadiusX];
					}
				}
				if (filterSum == 0) {
					pixels[y][x] = total;
				} else {
					pixels[y][x] = total / filterSum;
				}
			}
		}
		return { pixels, image.width, image.height };
	}

	void drawFilledCircle(ColorImage image, int x, int y, int radius, int* color) {
		for (int r = 0; r < radius; r++) {
			auto ring = getCirclePixels(x, y, r);
			for (const auto& pixel : ring) {
				int px = pixel.first;
				int py = pixel.second;
				if (inbounds(px, py, image.width, image.height)) {
					image.pixels[px][py] = color;
				}
			}
		}
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
			if (inbounds(x, y, image.width, image.height)) {
				image.pixels[y][x] = color;
			}
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
		GrayscaleImage newImage {
			new int*[magnitudes.height],
			magnitudes.width,
			magnitudes.height,
			1
		};

		for (int y = 0; y < magnitudes.height; y++) {
			newImage.pixels[y] = new int[magnitudes.width];
			for (int x = 0; x < magnitudes.width; x++) {
				int pixelValue = magnitudes.pixels[y][x];
				if (pixelValue > upperThreshold) {
					unvisited.insert({ x, y });
					newImage.pixels[y][x] = 2;
				} else if (pixelValue > lowerThreshold) {
					newImage.pixels[y][x] = 1;
				} else {
					newImage.pixels[y][x] = 0;
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
				if (x_ < 0 || x_ >= magnitudes.width) continue;
				for (int y_ = y - 1; y_ <= y + 1; y_++) {
					// Boundary check
					if (y_ < 0 || y_ >= magnitudes.height) continue;

					// If the current value is 1 (reached lower threshold, but not upper threshold),
					// then because it touches a strong edge with a value of 2, it is automatically
					// promoted. Then, we add it to the queue to update its neighbors as well.
					if (newImage.pixels[y_][x_] == 1) {
						newImage.pixels[y_][x_] = 2;
						unvisited.insert({x_, y_});
					}
				}
			}
		}

		for (int y = 0; y < newImage.height; y++) {
			for (int x = 0; x < newImage.width; x++) {
				// Convert 2 to 1, and 1 to 0.
				newImage.pixels[y][x] >>= 1;
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
		int width = xGradient.width;
		int height = xGradient.height;

		int **newPixels = new int*[height];
		for (int y = 0; y < height; y++) {
			newPixels[y] = new int[width];

			for (int x = 0; x < width; x++) {
				newPixels[y][x] = 0;

				double currentMagnitude = magnitudes.pixels[y][x];

				// angle is in the range [-pi, pi].
				double angleRadians = atan2(yGradient.pixels[y][x], xGradient.pixels[y][x]);
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
					if (currentMagnitude > magnitudes.pixels[y + dy][x + dx]) {
						if (inBounds(x - dx, y - dy, width, height)) {
							if (currentMagnitude > magnitudes.pixels[y - dy][x - dx]) {
								newPixels[y][x] = 1;
							}
						}
					}
				}
			}
		}

		return GrayscaleImage { newPixels, width, height, 1 };
	}

	GrayscaleImage combineSobel(GrayscaleImage first, GrayscaleImage second) {
		int** pixels = new int*[first.height];
		for (int y = 0; y < first.height; y++) {
			pixels[y] = new int[first.width];
			for (int x = 0; x < first.width; x++) {
				pixels[y][x] = (int) sqrt(findMagnitudeSquared(first.pixels[y][x], second.pixels[y][x]));
			}
		}

		return { pixels, first.width, first.height };
	}

	GrayscaleImage applyThreshold(GrayscaleImage image, double threshold) {
		int** pixels = new int*[image.height];
		for (int y = 0; y < image.height; y++) {
			pixels[y] = new int[image.width];
			for (int x = 0; x < image.width; x++) {
				pixels[y][x] = image.pixels[y][x] >= threshold ? 255 : 0;
			}
		}

		return { pixels, image.width, image.height };
	}

	GrayscaleImage combineImages(GrayscaleImage afterHysteresis, GrayscaleImage afterNonMaxSuppression) {
		GrayscaleImage newImage = { new int*[afterHysteresis.height], afterHysteresis.width, afterHysteresis.height, 1 };

		for (int y = 0; y < afterHysteresis.height; y++) {
			newImage.pixels[y] = new int[afterHysteresis.width];
			for (int x = 0; x < afterHysteresis.width; x++) {
				if (afterHysteresis.pixels[y][x] != 0 && afterNonMaxSuppression.pixels[y][x] != 0) {
					newImage.pixels[y][x] = 1;
				} else {
					newImage.pixels[y][x] = 0;
				}
			}
		}

		return newImage;
	}

	/**
	 * Returns a 2D array of radian angle measures with the same dimensions as the input images.
	 */
	double ** getEdgeAnglesFromGradients(GrayscaleImage xGradient, GrayscaleImage yGradient) {
		double ** measures = new double*[xGradient.height];
		for (int y = 0; y < xGradient.height; y++) {
			measures[y] = new double[xGradient.width];
			for (int x = 0; x < xGradient.width; x++) {
				measures[y][x] = atan2(yGradient.pixels[y][x], xGradient.pixels[y][x]);
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
		GrayscaleImage xGradient = convolve(grayscale, horizontalSobel);
		GrayscaleImage yGradient = convolve(grayscale, verticalSobel);
		
		GrayscaleImage magnitudes = combineSobel(xGradient, yGradient);

		GrayscaleImage afterHysteresis = hysteresis(magnitudes, lowerThreshold, upperThreshold);
		GrayscaleImage afterNonMaxSuppression = nonMaxSuppression(xGradient, yGradient, magnitudes);

		EdgeDetectionResult result;
		result.edges = combineImages(afterHysteresis, afterNonMaxSuppression);
		// tjcv::saveGrayscalePPM("image_h.ppm", afterHysteresis);
		result.angles = getEdgeAnglesFromGradients(xGradient, yGradient);
		result.xGradient = xGradient;
		result.yGradient = yGradient;
		
		return result;
	}
}

namespace lab6 {
	using tjcv::GrayscaleImage;

	/**
	 * Using Bresenham's algorithm, casts votes along a line. The rayWidth option allows you to cast lines that are wider than one pixel.
	 */
	void castVotesForOnePixel(int **votes, int x, int y, int width, int height, double angle, int rayWidth = 0) {
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
		int **votes = tjcv::make2DIntArray(edges.height, edges.width);

		for (int y = 0; y < edges.height; y++) {
			for (int x = 0; x < edges.width; x++) {
				if (edges.pixels[y][x]) {
					castVotesForOnePixel(votes, x, y, edges.width, edges.height, angles[y][x], rayWidth);
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
	 * countEdgesForCircle counts the number of pixels along the edge of a circle.
	 * The circle is specified by an x, y, and radius.
	 */
	int countEdgesForCircle(GrayscaleImage edges, int x, int y, int radius) {
		int count = 0;
		auto circlePixels = tjcv::getCirclePixels(x, y, radius);
		for (const auto& pixel : circlePixels) {
			int pixelX = pixel.first;
			int pixelY = pixel.second;
			if (tjcv::inbounds(pixelX, pixelY, edges.width, edges.height)) {
				if (edges.pixels[pixelY][pixelX] > 0) {
					count++;
				}
			}
		}

		return count;
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
		int *edgesByRadius = new int[maxRadius - minRadius];
		for (int radius = minRadius; radius < maxRadius; radius++) {
			edgesByRadius[radius - minRadius] = countEdgesForCircle(edges, x, y, radius);
		}
		
		std::vector<int> radii;

		// Iterate over each possible radius, using the given ringWidth.
		for (int radius = minRadius; radius < maxRadius; radius++) {
			// Keep track of how many edges are found along the circle
			int ringEdgeCount    = edgesByRadius[radius];
			// Keep track of how many edges could possibly be found along the circle
			int ringMaxEdgeCount = tjcv::countPixelsToDrawCircle(radius);

			// Ring radii must be in the interval requested
			int ringStartRadius  = __max(minRadius, radius - ringWidth);
			int ringEndRadius    = __min(radius + ringWidth, maxRadius - 1);

			for (int ringRadius = ringStartRadius; ringRadius <= ringEndRadius; ringRadius++) {
				ringMaxEdgeCount += tjcv::countPixelsToDrawCircle(ringRadius);
				ringEdgeCount    += edgesByRadius[ringRadius - minRadius];
			}

			double ratio = ((double) ringEdgeCount / ringMaxEdgeCount);
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

		GrayscaleImage image;
		image.width = width;
		image.height = height;
		image.pixels = new int*[height];
		image.max = maxIntensity;
		for (int y = 0; y < height; y++) {
			image.pixels[y] = new int[width];
			for (int x = 0; x < width; x++) {
				image.pixels[y][x] = votes[y][x];
			}
		}
		return image;
	}

	void part1() {
		dbg("Loading image\n");
		auto color = tjcv::loadColorPPM("image.ppm");
		auto grayscale = tjcv::convertToGrayscale(color);

		dbg("Detecting edges\n");
		auto detection = lab5::detectEdges(grayscale, 160, 180);
		tjcv::saveGrayscalePPM("imagef.ppm", detection.edges);

		dbg("Casting votes\n");
		auto votes = lab6::castVotes(detection.edges, detection.angles, 2);
		auto votesGraph = lab6::createVotesGraph(votes, detection.edges.width, detection.edges.height);
		tjcv::saveGrayscalePPM("imagev.ppm", votesGraph);

		dbg("Finding centers\n");
		auto centers = lab6::findCenters(votes, grayscale.width, grayscale.height, 200);

		dbg("Found " << centers.size() << " centers\n");

		dbg("Drawing centers\n");
		auto colorWithCenters = tjcv::cloneColorImage(color);
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
					if (tjcv::inbounds(px, py, color.width, color.height)) {
						colorWithCenters.pixels[py][px] = CENTER_COLOR;
					}
				}
			}
		}

		dbg("Saving centers\n");
		tjcv::saveColorPPM("imageCC.ppm", colorWithCenters);
	}

	void part2() {
		dbg("Loading image\n");
		auto color = tjcv::loadColorPPM("image.ppm");
		auto grayscale = tjcv::convertToGrayscale(color);

		dbg("Detecting edges\n");
		auto detection = lab5::detectEdges(grayscale, 160, 180);
		tjcv::saveGrayscalePPM("imagef.ppm", detection.edges);

		dbg("Casting votes\n");
		auto votes = lab6::castVotes(detection.edges, detection.angles, 2);
		auto votesGraph = lab6::createVotesGraph(votes, detection.edges.width, detection.edges.height);
		tjcv::saveGrayscalePPM("imagev.ppm", votesGraph);

		dbg("Finding centers\n");
		auto centers = lab6::findCenters(votes, grayscale.width, grayscale.height, 200);

		dbg("Found " << centers.size() << " centers\n");

		int* CIRCLE_COLOR = new int[3] { 0, 255, 0 };
		int* CENTER_COLOR = new int[3] { 255, 0, 0 };

		dbg("Finding radii\n");
		int foundRadiusCount = 0;
		for (int i = 0; i < centers.size(); i++) {
			const auto& center = centers.at(i);
			int x = center.first;
			int y = center.second;
			
			auto radii = lab6::findRadii(detection.edges, x, y, 50, 200, 2, 0.6);
			for (int radius : radii) {
				dbg("Circle: {x=" << x << ", y=" << y << ", r=" << radius << "}\n");
				tjcv::drawCircle(color, x, y, radius, CIRCLE_COLOR);
			}

			foundRadiusCount += radii.size();
		}

		dbg("Found " << foundRadiusCount << " radii\n");
		dbg("Saving\n");

		tjcv::saveColorPPM("imagecircles.ppm", color);
	}
}

int main() {
	lab6::part2();
}
