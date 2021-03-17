// l051.cpp
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
	} GrayscaleImage;

	typedef byte** Filter;
	
	/**
	 * hysteresis() takes an image and runs a double threshold on it. Then, any "weak" pixels
	 * touching a "strong" pixel are automatically promoted to "strong" pixels. This process is
	 * repeated until there are no more changes.
	 */
	GrayscaleImage hysteresis(GrayscaleImage image, int lowerThreshold, int upperThreshold) {
		// set of (x, y) pairs
		std::set<std::pair<int, int>> unvisited;
		GrayscaleImage newImage = GrayscaleImage { new byte*[image.height], image.width, image.height };

		for (int y = 0; y < image.height; y++) {
			newImage.pixels[y] = new byte[image.width];
			for (int x = 0; x < image.width; x++) {
				byte pixelValue = image.pixels[y][x];
				newImage.pixels[y][x] = (pixelValue > upperThreshold) + (pixelValue > lowerThreshold);
				if (pixelValue > upperThreshold) {
					unvisited.insert({ x, y });
				}
			}
		}

		while (unvisited.size()) {
			// Pop a location from the set
			auto next = *unvisited.begin();
			int x = next.first, y = next.second;
			unvisited.erase(unvisited.begin());

			// Check nearby locations. x_ and y_ represent points around the current pixel.
			for (int x_ = x - 1; x_ <= x + 1; x_++) {
				// Boundary check
				if (x_ < 0 || x_ >= image.width) continue;
				for (int y_ = y - 1; y_ <= y + 1; y_++) {
					// Boundary check
					if (y_ < 0 || y_ >= image.height) continue;

					// If the current value is 1 (reached lower threshold, but not upper threshold),
					// then because it touches a strong edge with a value of 2, it is automatically
					// promoted. Then, we add it to the queue to update its neighbors as well.
					if (image.pixels[y_][x_] == 1) {
						image.pixels[y_][x_] = 2;
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
	GrayscaleImage nonMaxSuppression(GrayscaleImage xGradient, GrayscaleImage yGradient) {
		int width = xGradient.width;
		int height = xGradient.height;

		double **magnitudesSquared = new double*[height];

		for (int y = 0; y < height; y++) {
			magnitudesSquared[y] = new double[width];
			for (int x = 0; x < width; x++) {
				magnitudesSquared[y][x] = findMagnitudeSquared(yGradient.pixels[y][x], xGradient.pixels[y][x]);
			}
		}

		byte **newPixels = new byte*[height];
		for (int y = 0; y < height; y++) {
			newPixels[y] = new byte[width];

			for (int x = 0; x < width; x++) {
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
				
				bool isMax = true;

				// Initialize the magnitudes, in case they are out of bounds.
				double currMagnitudeSquared = magnitudesSquared[y][x];

				if (inBounds(x + dx, y + dy, width, height)) {
					if (currMagnitudeSquared < magnitudesSquared[y + dy][x + dx]) {
						isMax = false;
					} else if (inBounds(x - dx, y - dy, width, height)) {
						if (currMagnitudeSquared < magnitudesSquared[y - dy][x - dx]) {
							isMax = false;
						}
					}
				}

				newPixels[y][x] = isMax ? ((byte) sqrt(currMagnitudeSquared)) : 0;
			}
		}

		return GrayscaleImage { newPixels, width, height };
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

		for (int y = 1; y < image.height - 1; y++) {
			for (int x = 1; x < image.width - 1; x++) {
				int total = 0;
				for (int relativeX = -1; relativeX <= 1; relativeX++) {
					for (int relativeY = -1; relativeY <= 1; relativeY++) {
						// std::cout << "Accessing pixel " << x + relativeX << ", " << y + relativeY << '\n';
						total += image.pixels[y + relativeY][x + relativeX] * filter[relativeY + 1][relativeX + 1];
					}
				}
				pixels[y][x] = total / 9;
			}
		}
		return { pixels, image.width, image.height };
	}

	GrayscaleImage combineSobel(GrayscaleImage first, GrayscaleImage second) {
		byte** pixels = new byte*[first.height];
		for (int y = 0; y < first.height; y++) {
			pixels[y] = new byte[first.width];
			for (int x = 0; x < first.width; x++) {
				pixels[y][x] = (int) sqrt(first.pixels[y][x] * first.pixels[y][x] + second.pixels[y][x] * second.pixels[y][x]);
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
				// std::cout << "(" << image.pixels[y][x][0] << ", " << image.pixels[y][x][1] << ", " << image.pixels[y][x][2] << ")\n";
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

	void part1() {
		ColorImage image = loadColorPPM("image.ppm");
		GrayscaleImage grayscale = convertToGrayscale(image);
		saveGrayscalePPM("imageg.ppm", grayscale);

		Filter horizontalSobel = new byte*[3] {
			new byte[3] {  1,  2,  1 },
			new byte[3] {  0,  0,  0 },
			new byte[3] { -1, -2, -1 }
		};
		Filter verticalSobel = new byte*[3] {
			new byte[3] {  1,  0, -1 },
			new byte[3] {  2,  0, -2 },
			new byte[3] {  1,  0, -1 }
		};

		GrayscaleImage horizontalFiltered = convolve(grayscale, horizontalSobel);
		GrayscaleImage verticalFiltered = convolve(grayscale, verticalSobel);
		GrayscaleImage combined = combineSobel(horizontalFiltered, verticalFiltered);
		GrayscaleImage thresholded = applyThreshold(combined, 45);

		saveGrayscalePPM("imagem.ppm", thresholded);

		// saveGrayscalePPM("imagehf.ppm", horizontalFiltered);
		// saveGrayscalePPM("imagevf.ppm", verticalFiltered);
	}
}


int main() {
	lab5::part1();	
}