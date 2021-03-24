// l053.cpp
// Michael Fatemi

#include <iostream>
#include <fstream>
#include <cmath>
#include <set>
#include <vector>

namespace lab5 {
	typedef struct {
		int*** pixels;
		int width, height;
	} ColorImage;

	typedef struct {
		int** pixels;
		int width, height;
	} GrayscaleImage;

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
		GrayscaleImage newImage = GrayscaleImage { new int*[magnitudes.height], magnitudes.width, magnitudes.height };

		for (int y = 0; y < magnitudes.height; y++) {
			newImage.pixels[y] = new int[magnitudes.width];
			for (int x = 0; x < magnitudes.width; x++) {
				int pixelValue = magnitudes.pixels[y][x];
				if (pixelValue > upperThreshold) {
					unvisited.insert({ x, y });
					newImage.pixels[y][x] = 255;
				} else if (pixelValue > lowerThreshold) {
					newImage.pixels[y][x] = 128;
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
					if (newImage.pixels[y_][x_] == 128) {
						newImage.pixels[y_][x_] = 255;
						unvisited.insert({x_, y_});
					}
				}
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

				if (currentMagnitude < threshold) {
					continue;
				}

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
								newPixels[y][x] = 255; // (int) sqrt(currMagnitudeSquared);
								// std::cout << currentMagnitude << " > {" << magnitudes.pixels[y + dy][x + dx] << ", " << magnitudes.pixels[y - dy][x - dx] << "}\n";
							}
						}
					}
				}
			}
		}

		return GrayscaleImage { newPixels, width, height };
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

	void saveGrayscalePPM(std::string filename, GrayscaleImage image) {
		std::ofstream handle(filename);
		// P2 [width] [height] [max intensity]
		handle << "P2 " << image.width << " " << image.height << " 255\n";
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
		GrayscaleImage newImage = { new int*[afterHysteresis.height], afterHysteresis.width, afterHysteresis.height };

		for (int y = 0; y < afterHysteresis.height; y++) {
			newImage.pixels[y] = new int[afterHysteresis.width];
			for (int x = 0; x < afterHysteresis.width; x++) {
				if (afterHysteresis.pixels[y][x] > 128 && afterNonMaxSuppression.pixels[y][x] != 0) {
					newImage.pixels[y][x] = 255;
				} else {
					newImage.pixels[y][x] = 0;
				}
			}
		}

		return newImage;
	}

	void part1() {
		ColorImage image = loadColorPPM("image.ppm");
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
		ColorImage image = loadColorPPM("image.ppm");
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

namespace graphicsutil {
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

	std::vector<std::pair<int, int>> projectCirclePixels(int centerX, int centerY, double radius) {
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
	
	int ** make2DIntArray(int dim1, int dim2) {
		int ** arr = new int*[dim1];
		for (int i = 0; i < dim1; i++) {
			arr[i] = new int[dim2];
		}

		return arr;
	}

	bool inbounds(int x, int y, int maxX, int maxY) {
		return (x >= 0 && x < maxX) && (y >= 0 && y < maxY);
	}
}

namespace lab6 {
	void castVotesForOnePixel(int ** votes, int x, int y, int width, int height, double angle) {
		using graphicsutil::inbounds;

		graphicsutil::BresenhamPixelIterator it(x, y, angle);
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

	int** castVotes(lab5::GrayscaleImage edges, double **angles) {
		int **votes = graphicsutil::make2DIntArray(edges.height, edges.width);

		for (int y = 0; y < edges.width; y++) {
			for (int x = 0; x < edges.height; x++) {
				if (edges.pixels[y][x]) {
					castVotesForOnePixel(votes, x, y, edges.width, edges.height, angles[y][x]);
				}
			}
		}

		return votes;
	}
}

int main() {
	
}
