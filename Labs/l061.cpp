// l053.cpp
// Michael Fatemi

#include <iostream>
#include <fstream>
#include <cmath>
#include <set>
#include <vector>

namespace tjcv {
	typedef struct {
		int*** pixels;
		int width, height;
	} ColorImage;

	typedef struct {
		int** pixels;
		int width, height;
		int max = 1;
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
					buildup += count * cos(angle) / sin(angle);
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
				handle << image.pixels[y][x] << " ";
			}
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

			pixels.push_back({ x + centerX, y + centerY });
			pixels.push_back({ x + centerX, -y + centerY });
			pixels.push_back({ -x + centerX, y + centerY });
			pixels.push_back({ -x + centerX, -y + centerY });

			pixels.push_back({ y + centerX, x + centerY });
			pixels.push_back({ y + centerX, -x + centerY });
			pixels.push_back({ -y + centerX, x + centerY });
			pixels.push_back({ -y + centerX, -x + centerY });

			y2_new -= 2 * x - 3;

			x += 1;
		}

		return pixels;
	}

	double ** make2DDoubleArray(int dim1, int dim2) {
		double ** arr = new double*[dim1];
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
}

namespace lab5 {
	using tjcv::GrayscaleImage;
	using tjcv::ColorImage;
	typedef struct {
		GrayscaleImage edges;
		double** angles;
	} EdgeDetectionResult;

	typedef int** Filter;

	Filter verticalSobel = nullptr, horizontalSobel = nullptr;
	
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
	GrayscaleImage nonMaxSuppression(GrayscaleImage xGradient, GrayscaleImage yGradient, GrayscaleImage magnitudes, int threshold) {
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

	GrayscaleImage convolve(GrayscaleImage image, Filter filter) {
		auto pixels = new int*[image.height];
		for (int y = 0; y < image.height; y++) {
			pixels[y] = new int[image.width];
			pixels[y][0] = 0;
			pixels[y][image.width - 1] = 0;
		}
		// Initialize the edges to 0
		for (int x = 0; x < image.width; x++) {
			pixels[0][x] = 0;
			pixels[image.height - 1][x] = 0;
		}

		for (int y = 1; y < image.height - 1; y++) {
			for (int x = 1; x < image.width - 1; x++) {
				int total = 0;
				for (int relativeX = -1; relativeX <= 1; relativeX++) {
					for (int relativeY = -1; relativeY <= 1; relativeY++) {
						total += image.pixels[y + relativeY][x + relativeX] * filter[relativeY + 1][relativeX + 1];
					}
				}
				pixels[y][x] = total / 9;
			}
		}
		return { pixels, image.width, image.height };
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

	void divideInPlace(GrayscaleImage image, int amount) {
		for (int y = 0; y < image.height; y++) {
			for (int x = 0; x < image.width; x++) {
				image.pixels[y][x] = (int) (image.pixels[y][x] / amount);
			}
		}
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

	void part1() {
		ColorImage image = tjcv::loadColorPPM("image.ppm");
		GrayscaleImage grayscale = convertToGrayscale(image);
		saveGrayscalePPM("imageg.ppm", grayscale);

		Filter verticalSobel = new int*[3] {
			new int[3] {	1,	2,	1 },
			new int[3] {	0,	0,	0 },
			new int[3] { -1, -2, -1 }
		};
		Filter horizontalSobel = new int*[3] {
			new int[3] {	1,	0, -1 },
			new int[3] {	2,	0, -2 },
			new int[3] {	1,	0, -1 }
		};

		GrayscaleImage verticalFiltered = convolve(grayscale, verticalSobel);
		GrayscaleImage horizontalFiltered = convolve(grayscale, horizontalSobel);
		GrayscaleImage combined = combineSobel(horizontalFiltered, verticalFiltered);
		GrayscaleImage thresholded = applyThreshold(combined, 45);

		saveGrayscalePPM("imagem.ppm", thresholded);
	}

	void part2() {
		ColorImage image = tjcv::loadColorPPM("image.ppm");
		GrayscaleImage grayscale = convertToGrayscale(image);

		Filter verticalSobel = new int*[3] {
			new int[3] {	1,	2,	1 },
			new int[3] {	0,	0,	0 },
			new int[3] { -1, -2, -1 }
		};
		Filter horizontalSobel = new int*[3] {
			new int[3] {	1,	0, -1 },
			new int[3] {	2,	0, -2 },
			new int[3] {	1,	0, -1 }
		};
		Filter gaussian = new int*[3] {
			new int[3] { 1, 2, 1 },
			new int[3] { 2, 4, 2 },
			new int[3] { 1, 2, 1 },
		};
	
		GrayscaleImage afterGaussian = convolve(grayscale, gaussian);
		// divideInPlace(afterGaussian, 16);
		GrayscaleImage xGradient = convolve(afterGaussian, horizontalSobel);
		GrayscaleImage yGradient = convolve(afterGaussian, verticalSobel);
		GrayscaleImage magnitudes = combineSobel(xGradient, yGradient);

		int lowerThreshold = 10;
		int upperThreshold = 30;

		GrayscaleImage afterHysteresis = hysteresis(magnitudes, lowerThreshold, upperThreshold);
		GrayscaleImage afterNonMaxSuppression = nonMaxSuppression(xGradient, yGradient, magnitudes, 10);
		GrayscaleImage finalResult = combineImages(afterHysteresis, afterNonMaxSuppression);

		// saveGrayscalePPM("image_magnitudes.ppm", magnitudes);
		saveGrayscalePPM("image1.ppm", afterNonMaxSuppression);
		saveGrayscalePPM("image2.ppm", afterHysteresis);
		saveGrayscalePPM("imagef.ppm", finalResult);
	}

	Filter getHorizontalSobel() {
		if (horizontalSobel == nullptr) {
			horizontalSobel = new int*[3] {
				new int[3] {	1,	0, -1 },
				new int[3] {	2,	0, -2 },
				new int[3] {	1,	0, -1 }
			};
		}

		return horizontalSobel;
	}

	Filter getVerticalSobel() {
		if (verticalSobel == nullptr) {
			verticalSobel = new int*[3] {
				new int[3] {	1,	2,	1 },
				new int[3] {	0,	0,	0 },
				new int[3] { -1, -2, -1 }
			};
		}

		return verticalSobel;
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
		Filter gaussian = new int*[3] {
			new int[3] { 1, 2, 1 },
			new int[3] { 2, 4, 2 },
			new int[3] { 1, 2, 1 },
		};
	
		GrayscaleImage afterGaussian = convolve(grayscale, gaussian);
		// divideInPlace(afterGaussian, 16);
		GrayscaleImage xGradient = convolve(afterGaussian, getHorizontalSobel());
		GrayscaleImage yGradient = convolve(afterGaussian, getVerticalSobel());
		GrayscaleImage magnitudes = combineSobel(xGradient, yGradient);

		GrayscaleImage afterHysteresis = hysteresis(magnitudes, lowerThreshold, upperThreshold);
		GrayscaleImage afterNonMaxSuppression = nonMaxSuppression(xGradient, yGradient, magnitudes, 10);

		EdgeDetectionResult result;
		result.edges = combineImages(afterHysteresis, afterNonMaxSuppression);
		result.angles = getEdgeAnglesFromGradients(xGradient, yGradient);
		
		return result;
	}
}

namespace lab6 {
	using tjcv::GrayscaleImage;

	/**
	 * Using Bresenham's algorithm, casts votes along a line. The rayWidth option allows you to cast lines that are wider than one pixel.
	 */
	void castVotesForOnePixel(int ** votes, int x, int y, int width, int height, double angle, int rayWidth = 0) {
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

	int** castVotes(GrayscaleImage edges, double **angles, int rayWidth = 0) {
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
	std::vector<int> findRadii(GrayscaleImage edges, int x, int y, int minRadius, int maxRadius, double minRatio) {
		std::vector<int> radii;
		for (int radius = minRadius; radius < maxRadius; radius++) {
			int edgesOnRadius = countEdgesForCircle(edges, x, y, radius);
			double maxEdgesOnRadius = radius * 2 * 3.141592653589;
			if ((edgesOnRadius / maxEdgesOnRadius) >= minRatio) {
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
		for (int y = 0; y < height; y++) {
			image.pixels[y] = new int[width];
			for (int x = 0; x < width; x++) {
				int intensity = (votes[y][x] * 255) / maxIntensity;
				image.pixels[y][x] = intensity < 255 ? intensity : 255;
			}
		}
		return image;
	}
}

int main() {
	auto color = tjcv::loadColorPPM("image.ppm");
	auto grayscale = tjcv::convertToGrayscale(color);
	std::cout << "Detecting edges\n";
	auto detection = lab5::detectEdges(grayscale, 10, 30);
	tjcv::saveGrayscalePPM("houghcircles_edges_output.ppm", detection.edges);
	std::cout << "Casting votes\n";
	auto votes = lab6::castVotes(detection.edges, detection.angles, 2);
	auto votesGraph = lab6::createVotesGraph(votes, detection.edges.width, detection.edges.height);
	tjcv::saveGrayscalePPM("houghcircles_votes_output.ppm", votesGraph);
	std::cout << "Finding centers\n";
	auto centers = lab6::findCenters(votes, grayscale.width, grayscale.height, 400);

	int* CIRCLE_COLOR = new int[3] { 0, 255, 0 };

	int foundRadiusCount = 0;
	std::cout << "Finding radii\n";
	std::cout << "Center count: " << centers.size() << '\n';
	for (int i = 0; i < centers.size(); i++) {
		const auto& center = centers.at(i);
		int x = center.first;
		int y = center.second;

		auto radii = lab6::findRadii(detection.edges, x, y, 50, 200, 0.3);
		for (int radius : radii) {
			std::cout << "Found radius " << radius << " for center " << (i + 1) << "/" << centers.size() << '\n';
			foundRadiusCount++;
			auto circle = tjcv::getCirclePixels(x, y, radius);
			for (auto circlePixel : circle) {
				int cpX = circlePixel.first;
				int cpY = circlePixel.second;
				if (tjcv::inbounds(cpX, cpY, color.width, color.height)) {
					color.pixels[cpY][cpX] = CIRCLE_COLOR;
				}
			}
		}
	}

	std::cout << "=== Summary ===\n";
	std::cout << "Found " << centers.size() << " centers,\n";
	std::cout << "Found " << foundRadiusCount << " radii\n";

	std::cout << "Saving\n";

	tjcv::saveColorPPM("houghcircles_color_output.ppm", color);
}
