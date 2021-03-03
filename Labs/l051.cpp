// l051.cpp
// Michael Fatemi

#include <iostream>
#include <fstream>
#include <cmath>

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