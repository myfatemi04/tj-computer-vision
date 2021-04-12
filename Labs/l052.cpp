// l052.cpp
// Michael Fatemi

#include <iostream>
#include <fstream>
#include <cmath>
#include <set>

namespace lab5 {
	typedef int byte;

	typedef struct {
		byte*** pixels;
		int width, height;
	} ColorImage;

	typedef struct {
		byte** pixels;
		int width, height;
		int max = 255;
	} GrayscaleImage;

	typedef byte** Filter;
	
	/**
	 * hysteresis() takes an image and runs a double threshold on it. Then, any "weak" edges
	 * touching a "strong" edge are automatically promoted to "strong" edges. This process is
	 * repeated until there are no more changes.
	 */
	GrayscaleImage hysteresis(GrayscaleImage magnitudes, int lowerThreshold, int upperThreshold) {
		// set of (x, y) pairs
		std::set<std::pair<int, int>> unvisited;
		GrayscaleImage newImage {
			new byte*[magnitudes.height],
			magnitudes.width,
			magnitudes.height,
			1
		};

		for (int y = 0; y < magnitudes.height; y++) {
			newImage.pixels[y] = new byte[magnitudes.width];
			for (int x = 0; x < magnitudes.width; x++) {
				byte pixelValue = magnitudes.pixels[y][x];
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

		byte **newPixels = new byte*[height];
		for (int y = 0; y < height; y++) {
			newPixels[y] = new byte[width];

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
		auto pixels = new byte*[image.height];
		for (int y = 0; y < image.height; y++) {
			pixels[y] = new byte[image.width];
			pixels[y][0] = 0;
			pixels[y][image.width - 1] = 0;
		}
		// Initialize the edges to 0
		for (int x = 0; x < image.width; x++) {
			pixels[0][x] = 0;
			pixels[image.height - 1][x] = 0;
		}

		int filterSum = 0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				filterSum += filter[i][j];
			}
		}

		for (int y = 1; y < image.height - 1; y++) {
			for (int x = 1; x < image.width - 1; x++) {
				int total = 0;
				for (int relativeX = -1; relativeX <= 1; relativeX++) {
					for (int relativeY = -1; relativeY <= 1; relativeY++) {
						total += image.pixels[y + relativeY][x + relativeX] * filter[relativeY + 1][relativeX + 1];
					}
				}
				pixels[y][x] = total / filterSum;
			}
		}
		return { pixels, image.width, image.height };
	}

	GrayscaleImage combineSobel(GrayscaleImage first, GrayscaleImage second) {
		byte** pixels = new byte*[first.height];
		for (int y = 0; y < first.height; y++) {
			pixels[y] = new byte[first.width];
			for (int x = 0; x < first.width; x++) {
				pixels[y][x] = (byte) sqrt(findMagnitudeSquared(first.pixels[y][x], second.pixels[y][x]));
			}
		}

		return { pixels, first.width, first.height };
	}

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
		auto pixels = new byte**[height];
		for (int y = 0; y < height; y++) {
			pixels[y] = new byte*[width];
			for (int x = 0; x < width; x++) {
				pixels[y][x] = new byte[3];
				handle >> pixels[y][x][0] >> pixels[y][x][1] >> pixels[y][x][2];
			}
		}
		handle.close();
		return { pixels, width, height };
	}

	GrayscaleImage convertToGrayscale(ColorImage image) {
		byte** pixels = new byte*[image.height];
		for (int y = 0; y < image.height; y++) {
			pixels[y] = new byte[image.width];
			for (int x = 0; x < image.width; x++) {
				pixels[y][x] = (image.pixels[y][x][0] + image.pixels[y][x][1] + image.pixels[y][x][2]) / 3;
			}
		}
		return { pixels, image.width, image.height };
	}

	GrayscaleImage applyThreshold(GrayscaleImage image, double threshold) {
		byte** pixels = new byte*[image.height];
		for (int y = 0; y < image.height; y++) {
			pixels[y] = new byte[image.width];
			for (int x = 0; x < image.width; x++) {
				pixels[y][x] = image.pixels[y][x] >= threshold ? 255 : 0;
			}
		}

		return { pixels, image.width, image.height };
	}

	void divideInPlace(GrayscaleImage image, byte amount) {
		for (int y = 0; y < image.height; y++) {
			for (int x = 0; x < image.width; x++) {
				image.pixels[y][x] = (byte) (image.pixels[y][x] / amount);
			}
		}
	}

	GrayscaleImage combineImages(GrayscaleImage afterHysteresis, GrayscaleImage afterNonMaxSuppression) {
		GrayscaleImage newImage = {
			new byte*[afterHysteresis.height], // pixels
			afterHysteresis.width,
			afterHysteresis.height,
			1 // max intensity
		};

		for (int y = 0; y < afterHysteresis.height; y++) {
			newImage.pixels[y] = new byte[afterHysteresis.width];
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
		ColorImage image = loadColorPPM("image.ppm");
		GrayscaleImage grayscale = convertToGrayscale(image);
		saveGrayscalePPM("imageg.ppm", grayscale);

		Filter verticalSobel = new byte*[3] {
			new byte[3] {  1,  2,  1 },
			new byte[3] {  0,  0,  0 },
			new byte[3] { -1, -2, -1 }
		};
		Filter horizontalSobel = new byte*[3] {
			new byte[3] {  1,  0, -1 },
			new byte[3] {  2,  0, -2 },
			new byte[3] {  1,  0, -1 }
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

		Filter verticalSobel = new byte*[3] {
			new byte[3] {  1,  2,  1 },
			new byte[3] {  0,  0,  0 },
			new byte[3] { -1, -2, -1 }
		};
		Filter horizontalSobel = new byte*[3] {
			new byte[3] {  1,  0, -1 },
			new byte[3] {  2,  0, -2 },
			new byte[3] {  1,  0, -1 }
		};
		Filter gaussian = new byte*[3] {
			new byte[3] { 1, 2, 1 },
			new byte[3] { 2, 4, 2 },
			new byte[3] { 1, 2, 1 },
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
}

int main() {
	lab5::part2();	
}
