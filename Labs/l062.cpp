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

int _divideRoundUp(int numerator, int denominator) {
	return (numerator + denominator - 1) / denominator;
}

double _safeAbs(double a) {
	return a < 0 ? -a : a;
}

int max(int a, int b) {
	return (a > b) ? a : b;
}

int min(int a, int b) {
	return (a < b) ? a : b;
}
namespace tjcv {
	/*
		This value is calculated by running the sobel operator along a perfect edge (0-->255).
	*/
	const double MAX_POSSIBLE_EDGE_GRADIENT = (2 * (127 + 127 * 2 + 127));

	enum ConvolveFillType {
		FILL_TYPE_BLACK,
		FILL_TYPE_WHITE,
		FILL_TYPE_ORIGINAL_PIXEL
	};

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
			int *get(int x, int y) const;

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
			int absoluteMax = 0;

		public:
			GrayscaleImage(int **pixels, int width, int height, int max = 255);
			GrayscaleImage(int width, int height, int max = 255);
	
			void save(std::string filename);
			
			void set(int x, int y, int value);
			int get(int x, int y) const;

			int getMax();
			int getWidth() const;
			int getHeight() const;
			std::pair<int, int> getSize();
			void setMax(int max);
			bool has(int x, int y) const;
			int getAbsoluteMax() const;

			ColorImage toColor();

			GrayscaleImage convolve(GrayscaleImage filter, int mode);

			GrayscaleImage applyGaussianFilter(int kernelSize, double standardDeviation, ConvolveFillType fillType);
			GrayscaleImage applySobel(int kernelSize, int mode);
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

	int *ColorImage::get(int x, int y) const {
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

	int GrayscaleImage::get(int x, int y) const {
		if (x < 0 || y < 0) return -1;
		if (x >= width || y >= height) return -1;
		return this->pixels[y][x];
	}

	int GrayscaleImage::getMax() { return this->max; }
	int GrayscaleImage::getWidth() const { return this->width; }
	int GrayscaleImage::getHeight() const { return this->height; }
	std::pair<int, int> GrayscaleImage::getSize() { return { this->width, this->height }; }

	void GrayscaleImage::setMax(int max) { this->max = max; }
	bool GrayscaleImage::has(int x, int y) const { return (x >= 0) && (y >= 0) && (x < this->width) && (y < this->height); }
	int GrayscaleImage::getAbsoluteMax() const {
		int max = 0, _val;
		for (int x = 0; x < getWidth(); x++) {
			for (int y = 0; y < getHeight(); y++) {
				if ((_val = get(x, y)) > max) {
					max = _val;
				}
			}
		}
		return max;
	}

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

	const int CONVOLVE_FILL_BLACK = 0;
	const int CONVOLVE_FILL_WHITE = 1;
	const int CONVOLVE_FILL_ORIGINAL_PIXEL_VALUE = 2;
	const double PI = 3.141592653589;
	const double E = 2.718281828;

	double _square(double a) {
		return a * a;
	}

	GrayscaleImage GrayscaleImage::convolve(GrayscaleImage filter, int fillMode = 0) {
		// dbg("Convolving an image around a filter " << filter.getWidth() << "x" << filter.getHeight() << '\n');

		GrayscaleImage convolved(width, height, max);
		int filterSum = 0;
		for (int i = 0; i < filter.getHeight(); i++) {
			for (int j = 0; j < filter.getWidth(); j++) {
				filterSum += filter.get(i, j);
			}
		}
		// dbg("Filter sum: " << filterSum << '\n');

		int filterRadiusX = filter.getWidth() >> 1;
		int filterRadiusY = filter.getHeight() >> 1;
		// dbg("Filter radius: " << filterRadiusX << ", " << filterRadiusY << '\n');

		// Initialize the edges to 0
		for (int y = 0; y < height; y++) {
			for (int i = 0; i < filterRadiusX; i++) {
				if (fillMode == 2) {
					convolved.set(i, y, get(i, y));
					convolved.set(width - i - 1, y, get(width - i - 1, y));
				} else {
					convolved.set(i, y, fillMode * getMax());
					convolved.set(width - i - 1, y, fillMode * getMax());
				}
			}
		}

		for (int x = 0; x < width; x++) {
			for (int i = 0; i < filterRadiusY; i++) {
				if (fillMode == 2) {
					convolved.set(x, i, get(x, i));
					convolved.set(x, height - i - 1, get(x, height - i - 1));
				} else {
					convolved.set(x, i, fillMode * getMax());
					convolved.set(x, height - i - 1, fillMode * getMax());
				}
			}
		}

		for (int y = filterRadiusY; y < height - filterRadiusY; y++) {
			for (int x = filterRadiusX; x < width - filterRadiusX; x++) {
				int total = 0;
				for (int relativeX = -filterRadiusX; relativeX <= filterRadiusX; relativeX++) {
					for (int relativeY = -filterRadiusY; relativeY <= filterRadiusY; relativeY++) {
						total += get(x + relativeX, y + relativeY) * filter.get(relativeX + filterRadiusX, relativeY + filterRadiusY);
					}
				}

				convolved.set(x, y, filterSum != 0 ? (total / filterSum) : total);
			}
		}
		
		return convolved;
	}

	GrayscaleImage GrayscaleImage::applyGaussianFilter(int kernelSize, double standardDeviation, ConvolveFillType fillType) {
		double **kernel = new double*[kernelSize];
		double prefix = 1 / (2 * PI * _square(standardDeviation));
		double kernelSum = 0;
		int k = (kernelSize - 1) / 2;
		for (int i = -k; i <= k; i++) {
			kernel[i + k] = new double[kernelSize];
			for (int j = -k; j <= k; j++) {
				double exponent = -(i * i + j * j) / (2 * _square(standardDeviation));
				kernel[i + k][j + k] = prefix * pow(E, exponent);
				kernelSum += kernel[i + k][j + k];
			}
		}
		// Normalize the kernel
		for (int i = -k; i <= k; i++) {
			for (int j = -k; j <= k; j++) {
				kernel[i + k][j + k] /= kernelSum;
			}
		}
		

		GrayscaleImage newImage (width, height, this->max);

		for (int x = 0; x < width; x++) {
			for (int y = 0; y < height; y++) {
				int valueIfEmpty = 0;
				if (fillType == FILL_TYPE_BLACK) {
					valueIfEmpty = 0;
				} else if (fillType == FILL_TYPE_ORIGINAL_PIXEL) {
					valueIfEmpty = get(x, y);
				}
				if (x < k || y < k) {
					newImage.set(x, y, valueIfEmpty);
				} else if (((x + k) >= width) || ((y + k) >= height)) {
					newImage.set(x, y, valueIfEmpty);
				} else {
					double totalValue = 0;
					for (int i = -k; i <= k; i++) {
						for (int j = -k; j <= k; j++) {
							double kernelValue = kernel[i + k][j + k];
							int pixelValue = get(x + i, y + j);
							totalValue += pixelValue * kernelValue;
						}
					}
					// dbg("totalValue: " << totalValue << ", originalValue: " << get(x, y) << '\n');
					newImage.set(x, y, totalValue);
				}
			}
		}
		
		return newImage;
	}

	const int SOBEL_TYPE_VERTICAL = 0;
	const int SOBEL_TYPE_HORIZONTAL = 1;

	GrayscaleImage GrayscaleImage::applySobel(int kernelSize, int type) {
		double **kernel = new double*[kernelSize];
		int k = (kernelSize - 1) / 2;
		for (int i = -k; i <= k; i++) {
			kernel[i + k] = new double[kernelSize];
			for (int j = -k; j <= k; j++) {
				if (i == 0 && j == 0) {
					kernel[i + k][j + k] = 0;
				} else {
					kernel[i + k][j + k] = (type == SOBEL_TYPE_HORIZONTAL ? i : j) / (i * i + j * j);
				}
			}
		}

		dbg("Created Sobel filter of type " << type << " for kernelSize " << kernelSize << '\n');

		GrayscaleImage newImage (width, height, this->max);

		for (int x = 0; x < width; x++) {
			for (int y = 0; y < height; y++) {
				if (y < k || x < k) {
					newImage.set(x, y, 0);
				} else if ((y + k >= height) || (x + k >= width)) {
					newImage.set(x, y, 0);
				} else {
					double totalValue = 0;
					for (int i = -k; i <= k; i++) {
						for (int j = -k; j < k; j++) {
							double kernelValue = kernel[i + k][j + k];
							int pixelValue = get(x + i, y + j);
							totalValue += pixelValue * kernelValue;
						}
					}
					newImage.set(x, y, totalValue);
				}
			}
		}

		dbg("Applied Sobel filter successfully\n");
	
		return newImage;
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
				
				iterateOverX = _safeAbs(cos(angle)) > _safeAbs(sin(angle));
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
		GrayscaleImage magnitudes;
		GrayscaleImage xGradient, yGradient;
	} EdgeDetectionResult;

	/**
	 * hysteresis() takes an image and runs a double threshold on it. Then, any "weak" edges
	 * touching a "strong" edge are automatically promoted to "strong" edges. This process is
	 * repeated until there are no more changes.
	 */
	GrayscaleImage hysteresis(GrayscaleImage magnitudes, int lowerThreshold, int upperThreshold, int gaussianFilterRadius) {
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

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				im.set(x, y, 0);
				double currentMagnitude = magnitudes.get(x, y);

				if (currentMagnitude == 0) {
					continue;
				}

				// angle is in the range [-pi, pi].
				double angleRadians = atan2(yGradient.get(x, y), xGradient.get(x, y));
				double angleDegrees = angleRadians / 3.1415926535897932 * 180;
				int nearestFortyFive = (int)(round(angleDegrees / 45) * 45);

				// Find the amount to go in the X and Y directions to follow the gradient

				int dx = 0, dy = 0;
				switch (nearestFortyFive) {
					// Exact opposites are equivalent
					case 45:
					case -135:
						dx = 1;
						dy = 1;
						break;
					case -45:
					case 135:
						dx = 1;
						dy = -1;
					case -90:
					case 90:
						dx = 0;
						dy = 1;
						break;
					case 0:
					case 180:
					case -180:
						dx = 1;
						dy = 0;
						break;
				}

				int nextMagnitude = magnitudes.get(x + dx, y + dy);
				int prevMagnitude = magnitudes.get(x - dx, y - dy);

				im.set(x, y, (currentMagnitude >= nextMagnitude) && (currentMagnitude >= prevMagnitude));
			}
		}
		
		return im;
	}

	GrayscaleImage combineSobel(GrayscaleImage first, GrayscaleImage second) {
		GrayscaleImage combined(first.getWidth(), first.getHeight(), first.getMax());
		int max = 1;
		for (int y = 0; y < first.getHeight(); y++) {
			for (int x = 0; x < first.getWidth(); x++) {
				int value = (int) sqrt(findMagnitudeSquared(first.get(x, y), second.get(x, y)));
				combined.set(x, y, value);
				if (value > max) {
					max = value;
				}
			}
		}
		combined.setMax(max);
		
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
	EdgeDetectionResult detectEdges(GrayscaleImage grayscale, int lowerThreshold, int upperThreshold, int gaussianFilterRadius) {

		GrayscaleImage afterGaussian = grayscale.applyGaussianFilter(2 * gaussianFilterRadius + 1, 3, tjcv::FILL_TYPE_ORIGINAL_PIXEL);
		afterGaussian.save("imageg.ppm");
		GrayscaleImage xGradient = afterGaussian.applySobel(9, tjcv::SOBEL_TYPE_HORIZONTAL);
		GrayscaleImage yGradient = afterGaussian.applySobel(9, tjcv::SOBEL_TYPE_VERTICAL);
		
		GrayscaleImage magnitudes = combineSobel(xGradient, yGradient);

		GrayscaleImage afterHysteresis = hysteresis(magnitudes, lowerThreshold, upperThreshold, gaussianFilterRadius);
		GrayscaleImage afterNonMaxSuppression = nonMaxSuppression(xGradient, yGradient, magnitudes);
		afterHysteresis.save("imageh.ppm");
		afterNonMaxSuppression.save("imagenms.ppm");

		return EdgeDetectionResult {
			combineImages(afterHysteresis, afterNonMaxSuppression),
			getEdgeAnglesFromGradients(xGradient, yGradient),
			magnitudes,
			xGradient,
			yGradient
		};
	}
}

namespace lab6 {
	using tjcv::GrayscaleImage;

	typedef struct {
		int x, y, value;
	} _Maximum2D;

	typedef struct {
		int radius;
		double score;
	} RadiusResult;

	typedef struct {
		int x, y;
		RadiusResult radius;
	} CircleResult;

	/** Required functions to make sets of these types possible **/

	bool operator<(const RadiusResult& a, const RadiusResult& b) {
		return (a.radius < b.radius) || (a.score < b.score);
	}

	bool operator<(const CircleResult& a, const CircleResult& b) {
		return (a.x < b.x) || (a.y < b.y) || (a.radius < b.radius);
	}


	_Maximum2D __findRowMaximum(GrayscaleImage values, int width, int __x, int y) {
		int maximumValue = -1;
		int maximumValueX = -1;
		int maximumValueY = -1;
		for (int x = __x; x < __x + width && x < values.getWidth(); x++) {
			int value = values.get(x, y);
			if (value > maximumValue) {
				maximumValue = value;
				maximumValueX = x;
				maximumValueY = y;
			}
		}
		return _Maximum2D {maximumValueX, maximumValueY, maximumValue};
	}

	_Maximum2D __getMaximumFromCache(_Maximum2D *cache, int cacheSize) {
		_Maximum2D maximum {-1, -1, -1};
		for (int i = 0; i < cacheSize; i++) {
			auto current = cache[i];
			if (current.value > maximum.value) {
				maximum.value = current.value;
				maximum.x = current.x;
				maximum.y = current.y;
			}
		}
		return maximum;
	}

	std::set<std::pair<int, int>> findLocalMaximaWithSlidingSquare(GrayscaleImage values, int squareSize) {
		std::set<std::pair<int, int>> localMaxima;
	
		/*
		This caches the maximum value of each row of the current square.
		When the square moves down, the next value replaces the value one squareWidth
		before it in the cache. So, rowMaximumCache[y % squareSize] will be consistent.
		*/
		_Maximum2D *rowMaximumCache = new _Maximum2D[squareSize];

		for (int x = 0; x + squareSize < values.getWidth(); x++) {
			// Generate row maximum cache
			for (int y = 0; y < squareSize && y < values.getWidth(); y++) {
				auto rowMaximum = __findRowMaximum(values, squareSize, x, y);
				rowMaximumCache[y] = rowMaximum;
			}

			for (int nextRowY = squareSize; nextRowY < values.getHeight(); nextRowY++) {
				_Maximum2D rowMaximum = __findRowMaximum(values, squareSize, x, nextRowY);
				int rowCacheIndex = nextRowY % squareSize;
				rowMaximumCache[rowCacheIndex] = rowMaximum;
				// Now, we want to avoid the situation where the maximum value of the square keeps
				// increasing as the square goes down the image.
				// So, we only check the center row of the square.
				bool centerRowIsMax = true;
				int thisRowY = nextRowY - squareSize;
				int centerRowY = thisRowY + squareSize / 2;
				int centerColX = x + squareSize / 2;
				_Maximum2D centerRowMaximum = rowMaximumCache[centerRowY % squareSize];
				
				bool centerRowLocalMaximumIsAlsoCenterColumn = centerColX == centerRowMaximum.x;
				if (!centerRowLocalMaximumIsAlsoCenterColumn) {
					continue;
				}

				for (int testRowY = thisRowY; testRowY <= nextRowY; testRowY++) {
					int testRowCacheIndex = testRowY % squareSize;
					if (rowMaximumCache[testRowCacheIndex].value > centerRowMaximum.value) {
						centerRowIsMax = false;
					}
				}
				if (centerRowIsMax) {
					localMaxima.insert({centerRowMaximum.x, centerRowMaximum.y});
				}
				// _Maximum2D squareMaximum = __getMaximumFromCache(rowMaximumCache, squareSize);
				// localMaxima.insert({squareMaximum.x, squareMaximum.y});
			}
		}
		return localMaxima;
	}

	/**
	 * Using Bresenham's algorithm, casts votes along a line. The rayWidth option allows you to cast lines that are wider than one pixel.
	 */
	void castVotesForOnePixel(
		GrayscaleImage votes,
		int x,
		int y,
		double angle,
		int maxRadius) {
		using tjcv::inbounds;

		tjcv::BresenhamPixelIterator it(x, y, angle);

		int i, cx, cy;

		i = 0;
		while (cx = it.getX(), cy = it.getY(), votes.has(cx, cy) && (i < maxRadius || maxRadius == -1)) {
			votes.set(cx, cy, votes.get(cx, cy) + 1);
			it.step(1);
			i++;
		}

		it.reset();

		i = 0;
		while (cx = it.getX(), cy = it.getY(), votes.has(cx, cy) && (i < maxRadius || maxRadius == -1)) {
			votes.set(cx, cy, votes.get(cx, cy) + 1);
			it.step(-1);
			i++;
		}
	}

	GrayscaleImage castVotes(GrayscaleImage edges, double **angles, int maxRadius) {
		GrayscaleImage votes(edges.getWidth(), edges.getHeight());

		for (int y = 0; y < edges.getHeight(); y++) {
			for (int x = 0; x < edges.getWidth(); x++) {
				if (edges.get(x, y)) {
					castVotesForOnePixel(votes, x, y, angles[y][x], maxRadius);
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
	std::set<std::pair<int, int>> findCenters(GrayscaleImage votes, int localMaximumSlidingSquareSize, int threshold) {
		std::set<std::pair<int, int>> localMaxima = findLocalMaximaWithSlidingSquare(votes, localMaximumSlidingSquareSize);
		std::set<std::pair<int, int>> centers;
		for (const auto& localMaximum : localMaxima) {
			int x = localMaximum.first;
			int y = localMaximum.second;
			if (threshold != -1 && votes.get(x, y) > threshold) {
				centers.insert({x, y});
			}
		}

		return centers;
	}

	/**
	 * scoreCircleCandidate scores the likelihood of a radius being a circle.
	 * The way it does this is by iterating over the pixels that would compose a circle.
	 * The magnitude of the pixel's value multiplied by how well it's angled towards the center
	 * are accumulated for each pixel to produce a total score.
	 */
	double scoreCircleCandidate(
		GrayscaleImage magnitudes,
		double **angles,
		int x,
		int y,
		int radius,
		int maxPossibleEdgeGradient) {
		double score = 0;
		auto pixels = tjcv::getCirclePixels(x, y, radius);

		for (const auto& pixel : pixels) {
			int pixelX = pixel.first;
			int pixelY = pixel.second;
			
			// Calculate the dot product between the line along the edge gradient
			// and the ray coming from the center of the circle.

			// Angle
			double intendedAngle = atan2((double) (pixelY - y), (double) (pixelX - x));
			double actualAngle = angles[y][x];
			double theta = actualAngle - intendedAngle;
			double angleMatch = _safeAbs(cos(theta));

			// Magnitude
			int magnitude = magnitudes.get(pixelX, pixelY);
			double magnitudeScaled = (double) magnitude / maxPossibleEdgeGradient;

			if (magnitudeScaled > 0.1) {
				magnitudeScaled = 1;
			}

			// The magnitude from the center of the circle is scaled to be 1.
			double dotProduct = magnitudeScaled * angleMatch;

			// For some reason, the absolute value function outputs '0'
			score += dotProduct < 0 ? -dotProduct : dotProduct;
		}

		return score;
	}

	/**
	 * findRadii returns the radii that contain a certain number of edge pixels along their circumference. The radius is checked on the interval [minRadius, maxRadius). minScore specifies the minimum ratio of empty edges to fillled edges.
	 */
	std::vector<RadiusResult> findRadii(
		GrayscaleImage magnitudes,
		double **angles,
		int x,
		int y,
		int minRadius,
		int maxRadius,
		int ringWidth,
		double minScore,
		int maxPossibleEdgeGradient) {
		// Count the edges on each radius
		double *radiusScores = new double[maxRadius - minRadius];
		for (int radius = minRadius; radius < maxRadius; radius++) {
			radiusScores[radius - minRadius] = scoreCircleCandidate(magnitudes, angles, x, y, radius, maxPossibleEdgeGradient);
		}

		std::vector<RadiusResult> radii;

		// Iterate over each possible radius, using the given ringWidth.
		for (int radius = minRadius; radius < maxRadius; radius++) {
			// Keep track of how many edges are found along the circle
			int ringCumulativeScore    = radiusScores[radius - minRadius];
			// Keep track of how many edges could possibly be found along the circle
			int ringMaximumCumulativeScore = tjcv::countPixelsToDrawCircle(radius);

			// Ring radii must be in the interval requested
			// int ringStartRadius = radius;
			// int ringEndRadius   = min(radius + 1, maxRadius - 1);
			int ringStartRadius  = max(minRadius, radius - ringWidth);
			int ringEndRadius    = min(radius + ringWidth, maxRadius - 1);

			for (int ringRadius = ringStartRadius; ringRadius <= ringEndRadius; ringRadius++) {
				// ringMaxEdgeCount += tjcv::countPixelsToDrawCircle(ringRadius);
				ringCumulativeScore    += radiusScores[ringRadius - minRadius];
			}

			double ratio = (double) ringCumulativeScore / ringMaximumCumulativeScore;
			if (ratio > minScore) {
				radii.push_back({radius, ratio});
			} else if (ratio > minScore / 2) {
				// dbg("Found almost circle with ratio " << ratio << '\n');
			}
		}
		return radii;
	}

	/**
	 * When there are several circles in roughly the same area, we use deduplication.
	 * The way this works is by dividing the image into a grid. Then, we iterate through
	 * the circles, keeping track of which grid squares already have a circle. If a grid
	 * square already has a circle, we won't add it to the deduplicated result.
	 */
	std::set<CircleResult> dedupe(const std::set<CircleResult>& circles, int gridSquareWidth, int gridSquareHeight, int imageWidth, int imageHeight) {
		std::set<CircleResult> deduplicated;
		int gridSubsquareWidth = gridSquareWidth / 3;
		int gridSubsquareHeight = gridSquareHeight / 3;
		int horizontalGridSquareCount = _divideRoundUp(imageWidth, gridSubsquareWidth);
		int verticalGridSquareCount   = _divideRoundUp(imageHeight, gridSubsquareHeight);

		const CircleResult* **occupiedGridSquares = new const CircleResult* *[verticalGridSquareCount];
		for (int y = 0; y < verticalGridSquareCount; y++) {
			occupiedGridSquares[y] = new const CircleResult* [horizontalGridSquareCount];
			for (int x = 0; x < horizontalGridSquareCount; x++) {
				occupiedGridSquares[y][x] = nullptr;
			}
		}

		for (auto& circle : circles) {
			int x = circle.x;
			int y = circle.y;
			int gridSquareX = x / gridSubsquareWidth;
			int gridSquareY = y / gridSubsquareHeight;

			int maxScoreInSurroundingGridSquare = -1;
			int maxScoreSurroundingGridSquareX = -1;
			int maxScoreSurroundingGridSquareY = -1;

			// Check surrounding grid squares
			for (int i = gridSquareY - 1; i <= gridSquareY + 1; i++) {
				if (i < 0 || i >= verticalGridSquareCount) continue;
				for (int j = gridSquareX - 1; j <= gridSquareX + 1; j++) {
					if (j < 0 || j >= horizontalGridSquareCount) continue;

					const CircleResult *gridSquareResult = occupiedGridSquares[i][j];
					if (gridSquareResult == nullptr) {
						continue;
					}

					int scoreInGridSquare = gridSquareResult->radius.score;
					if (scoreInGridSquare > maxScoreInSurroundingGridSquare) {
						maxScoreInSurroundingGridSquare = scoreInGridSquare;
						maxScoreSurroundingGridSquareY = i;
						maxScoreSurroundingGridSquareX = j;
					}
				}
			}

			if (circle.radius.score > maxScoreInSurroundingGridSquare) {
				if (maxScoreSurroundingGridSquareX >= 0 && maxScoreSurroundingGridSquareY >= 0) {
					auto removedCircle = occupiedGridSquares[maxScoreSurroundingGridSquareY][maxScoreSurroundingGridSquareX];
					occupiedGridSquares[maxScoreSurroundingGridSquareY][maxScoreSurroundingGridSquareX] = nullptr;
					dbg("Removing a surrounding circle because this one had a higher score.\n");
					dbg("The other circle had x = " << removedCircle->x << ", y = " << removedCircle->y << '\n');
					dbg("This circle had x = " << x << ", y = " << y << '\n');
				}
				occupiedGridSquares[gridSquareY][gridSquareX] = new CircleResult { x, y, circle.radius };
			}
		}

		for (int gridSquareY = 0; gridSquareY < verticalGridSquareCount; gridSquareY++) {
			for (int gridSquareX = 0; gridSquareX < horizontalGridSquareCount; gridSquareX++) {
				const CircleResult *gridSquareResult = occupiedGridSquares[gridSquareY][gridSquareX];
				if (gridSquareResult != nullptr) {
					deduplicated.insert(*gridSquareResult);
				}
			}
		}

		return deduplicated;
	}

	std::set<std::pair<int, int>> dedupeCenters(const std::set<std::pair<int, int>>& centers, int gridSquareWidth, int gridSquareHeight, int imageWidth, int imageHeight) {
		std::set<std::pair<int, int>> deduplicated;
		int horizontalGridSquareCount = _divideRoundUp(imageWidth, gridSquareWidth);
		int verticalGridSquareCount   = _divideRoundUp(imageHeight, gridSquareHeight);

		bool **occupiedGridSquares = new bool*[verticalGridSquareCount];
		for (int y = 0; y < verticalGridSquareCount; y++) {
			occupiedGridSquares[y] = new bool[horizontalGridSquareCount];
			for (int x = 0; x < horizontalGridSquareCount; x++) {
				occupiedGridSquares[y][x] = false;
			}
		}

		for (const auto& center : centers) {
			int x = center.first;
			int y = center.second;
			int gridSquareX = x / gridSquareWidth;
			int gridSquareY = y / gridSquareHeight;

			if (!occupiedGridSquares[gridSquareY][gridSquareX]) {
				occupiedGridSquares[gridSquareY][gridSquareX] = true;
				deduplicated.insert(center);
			}
		}

		return deduplicated;
	}

	typedef struct {
		int pennies, nickels, dimes, quarters;
	} MoneyCountingResult;

	enum CoinType {
		COIN_TYPE_PENNY = 0,
		COIN_TYPE_NICKEL = 1,
		COIN_TYPE_DIME = 2,
		COIN_TYPE_QUARTER = 3,
		COIN_TYPE_HALF_DOLLAR = 4,
		COIN_TYPE_SILVER_DOLLAR = 5
	};

	double absoluteError(double a, double b) {
		return _safeAbs(a - b);
	}

	const double PENNY_RADIUS_MM = 19.05;

	const double* COIN_SIZE_RELATIVE_TO_PENNY = new double[6] {
		1, // penny
		21.21 / PENNY_RADIUS_MM, // nickel
		17.91 / PENNY_RADIUS_MM, // dime
		24.26 / PENNY_RADIUS_MM, // quarter
		30.61 / PENNY_RADIUS_MM, // half dollar
		38.1 / PENNY_RADIUS_MM   // silver dollar
	};

	const CoinType* coinTypes = new CoinType[6] {
		COIN_TYPE_PENNY,
		COIN_TYPE_NICKEL,
		COIN_TYPE_DIME,
		COIN_TYPE_QUARTER,
		COIN_TYPE_HALF_DOLLAR,
		COIN_TYPE_SILVER_DOLLAR
	};

	CoinType estimateCoinType(double unknownCoinRadius, double pennyRadius) {
		double radiusRelativeToPenny = unknownCoinRadius / pennyRadius;

		CoinType leastErrorCoinType = COIN_TYPE_PENNY;
		double leastErrorCoinTypeError = absoluteError(1, unknownCoinRadius);
		
		for (int i = 0; i < 6; i++) {
			double error = absoluteError(radiusRelativeToPenny, COIN_SIZE_RELATIVE_TO_PENNY[i]);
			if (error < leastErrorCoinTypeError) {
				leastErrorCoinTypeError = error;
				leastErrorCoinType = coinTypes[i];
			}
		}

		return leastErrorCoinType;
	}

	double findPennyRadius(const std::set<CircleResult>& circles, const tjcv::ColorImage& reference) {
		dbg("Finding penny radius\n");

		double pennyRadiiCumulativeTotal = 0;
		int pennyCount = 0;
		for (const auto& circle : circles) {
			auto projectedPixels = tjcv::getCirclePixels(circle.x, circle.y, circle.radius.radius / 2);
			int redPixelCount = 0;
			int totalPixelCount = projectedPixels.size();
			for (const auto& location : projectedPixels) {
				int x = location.first;
				int y = location.second;
				auto color = reference.get(x, y);
				int red = color[0];
				int green = color[1];
				int blue = color[2];
				double averageColor = (red + green + blue) / 3.0;
				if (red > (averageColor * 1.2)) {
					redPixelCount++;
				}
			}
			if (redPixelCount > totalPixelCount * 0.3) {
				pennyCount++;
				pennyRadiiCumulativeTotal += circle.radius.radius;
			}
		}
		if (pennyCount == 0) {
			return -1;
		} else {
			return pennyRadiiCumulativeTotal / pennyCount;
		}
	}

	void part2() {
		using tjcv::ColorImage;

		dbg("Loading image\n");
		auto colorImage = ColorImage::fromPPM("image.ppm");
		auto grayscaleImage = colorImage.toGrayscale();

		const int MIN_RADIUS = 60;
		const int MAX_RADIUS = 200;
		const int VOTE_LENGTH = -1;//MAX_RADIUS;
		const double SCORE_THRESHOLD = 1.5;
		const int RING_WIDTH = 2;
		
		// 2 * EDGES_GAUSSIAN_KERNEL_RADIUS + 1 = Kernel Size
		const int EDGES_GAUSSIAN_KERNEL_RADIUS = 14;
		const int MAX_POSSIBLE_EDGE_GRADIENT = tjcv::MAX_POSSIBLE_EDGE_GRADIENT / EDGES_GAUSSIAN_KERNEL_RADIUS;
		const int EDGE_LOWER_THRESHOLD = MAX_POSSIBLE_EDGE_GRADIENT * 0.1;
		const int EDGE_UPPER_THRESHOLD = MAX_POSSIBLE_EDGE_GRADIENT * 0.2;

		dbg("Edge lower threshold: " << EDGE_LOWER_THRESHOLD << '\n');
		dbg("Edge upper threshold: " << EDGE_UPPER_THRESHOLD << '\n');

		const int CENTER_LOCAL_MAXIMUM_SQUARE_SIZE = MIN_RADIUS * 2;
		const int CIRCLE_DEDUPE_SQUARE_SIZE = 25;
		const int CENTER_VOTES_THRESHOLD = 10;//3 * SCALE_FACTOR;
		const int VOTES_GAUSSIAN_KERNEL_SIZE = 9;
		const int VOTES_GAUSSIAN_STDDEV = 3;

		dbg("Detecting edges\n");
		auto detection = lab5::detectEdges(grayscaleImage, EDGE_LOWER_THRESHOLD, EDGE_UPPER_THRESHOLD, EDGES_GAUSSIAN_KERNEL_RADIUS);
		detection.edges.save("imagef.ppm");
		detection.magnitudes.save("imagemagnitudes.ppm");

		dbg("Casting votes\n");
		auto votes = lab6::castVotes(detection.edges, detection.angles, VOTE_LENGTH);
		votes = votes.applyGaussianFilter(VOTES_GAUSSIAN_KERNEL_SIZE, VOTES_GAUSSIAN_STDDEV, tjcv::FILL_TYPE_BLACK);
		votes.setMax(votes.getAbsoluteMax());
		votes.save("imagev.ppm");

		dbg("Finding centers\n");
		auto centers = lab6::findCenters(votes, CENTER_LOCAL_MAXIMUM_SQUARE_SIZE, CENTER_VOTES_THRESHOLD);
		// centers = dedupeCenters(centers, CENTER_DEDUPE_SQUARE_SIZE, CENTER_DEDUPE_SQUARE_SIZE, colorImage.getWidth(), colorImage.getHeight());

		dbg("Found " << centers.size() << " centers\n");

		int* CENTER_COLOR = new int[3] { 255, 0, 0 };

		auto imageWithCenters = colorImage.clone();

		std::set<CircleResult> tentativeCircles, deduplicatedCircles;

		double highestRadiusScore = 0;

		dbg("Finding radii\n");
		for (const auto& center : centers) {
			int x = center.first;
			int y = center.second;
			tjcv::drawFilledCircle(imageWithCenters, x, y, 5, CENTER_COLOR);
			
			auto radii = lab6::findRadii(detection.magnitudes, detection.angles, x, y, MIN_RADIUS, MAX_RADIUS, RING_WIDTH, SCORE_THRESHOLD, MAX_POSSIBLE_EDGE_GRADIENT);
			int bestRadius = -1;
			double bestRadiusScore = -1;
			for (const auto& radiusResult : radii) {
				if (radiusResult.score > bestRadiusScore) {
					bestRadiusScore = radiusResult.score;
					bestRadius = radiusResult.radius;
				}
			}
			if (bestRadiusScore > 0) {
				tentativeCircles.insert({x, y, RadiusResult {bestRadius, bestRadiusScore}});
				if (bestRadiusScore > highestRadiusScore) {
					highestRadiusScore = bestRadiusScore;
				}
			}
		}

		dbg("Deduplicating\n");
		deduplicatedCircles = tentativeCircles;
		deduplicatedCircles = dedupe(tentativeCircles, CIRCLE_DEDUPE_SQUARE_SIZE, CIRCLE_DEDUPE_SQUARE_SIZE, colorImage.getWidth(), colorImage.getHeight());

		dbg("Deduplicated finished successfully\n");

		int** CIRCLE_COLORS = new int*[6];
		CIRCLE_COLORS[COIN_TYPE_PENNY] = new int[3] { 255, 0, 0 };
		CIRCLE_COLORS[COIN_TYPE_NICKEL] = new int[3] { 255, 0, 255 };
		CIRCLE_COLORS[COIN_TYPE_DIME] = new int[3] { 0, 0, 255 };
		CIRCLE_COLORS[COIN_TYPE_QUARTER] = new int[3] { 0, 255, 0 };
		CIRCLE_COLORS[COIN_TYPE_SILVER_DOLLAR] = new int[3] { 255, 255, 0 };
		
		int *coinValuesByType = new int[6];
		coinValuesByType[COIN_TYPE_PENNY] = 1;
		coinValuesByType[COIN_TYPE_NICKEL] = 5;
		coinValuesByType[COIN_TYPE_DIME] = 10;
		coinValuesByType[COIN_TYPE_QUARTER] = 25;
		coinValuesByType[COIN_TYPE_SILVER_DOLLAR] = 100;

		double pennyRadius = findPennyRadius(deduplicatedCircles, colorImage);

		dbg("Found penny radius to be " << pennyRadius << " pixels\n");

		int *coinCountsByType = new int[6];
		for (int i = 0; i < 6; i++) coinCountsByType[i] = 0;

		int totalValueInCents = 0;

		for (const auto& circle : deduplicatedCircles) {
			double coinRadius = circle.radius.radius;

			CoinType estimatedCoinType = estimateCoinType(coinRadius, pennyRadius);
			int* circleColor = CIRCLE_COLORS[estimatedCoinType];

			coinCountsByType[estimatedCoinType]++;
			totalValueInCents += coinValuesByType[estimatedCoinType];

			tjcv::drawCircle(colorImage, circle.x, circle.y, circle.radius.radius, circleColor);
			dbg("Coin: {x=" << circle.x << ", y=" << circle.y << ", r=" << coinRadius << ", score=" << circle.radius.score << "}\n");
		}

		imageWithCenters.save("imageCC.ppm");

		dbg("Found " << tentativeCircles.size() << " radii before deduplication\n");
		dbg("Found " << deduplicatedCircles.size() << " radii after deduplication\n");
		dbg("Saving\n");

		colorImage.save("coins.ppm");

		int dollarCount = totalValueInCents / 100;
		int centCount = totalValueInCents % 100;

		std::cout << "Result:\n";
		std::cout << coinCountsByType[COIN_TYPE_PENNY] << " pennies\n";
		std::cout << coinCountsByType[COIN_TYPE_NICKEL] << " nickels\n";
		std::cout << coinCountsByType[COIN_TYPE_DIME] << " dimes\n";
		std::cout << coinCountsByType[COIN_TYPE_QUARTER] << " quarters\n";
		std::cout << coinCountsByType[COIN_TYPE_SILVER_DOLLAR] << " silver dollars\n";
		std::cout << "For a grand total of $" << dollarCount << "." << (centCount < 10 ? "0" : "") << centCount << "!\n";
	}
}

int main() {
	lab6::part2();
}
