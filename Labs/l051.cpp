#include <iostream>
#include <fstream>

namespace lab5 {
	typedef unsigned char byte;

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
		// Initialize the edges to 0
		for (int x = 0; x < image.width; x++) {
			pixels[0][x] = 0;
			pixels[image.height - 1][x] = 0;
		}
		for (int y = 0; y < image.height; y++) {
			pixels[y][0] = 0;
			pixels[y][image.width - 1] = 0;
		}

		for (int y = 1; y < image.height - 1; y++) {
			pixels[y] = new byte[image.width];
			for (int x = 1; x < image.width - 1; x++) {
				int total = 0;
				for (int relativeX = -1; relativeX <= 1; relativeX++) {
					for (int relativeY = -1; relativeY <= 1; relativeY++) {
						total += image.pixels[y + relativeY][x + relativeX];
					}
				}
				byte value = total / 9;
				pixels[y][x] = value;
			}
		}
		return { pixels, image.width, image.height };
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
		int width, height;
		handle >> width >> height;
		int _max;
		handle >> _max;
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

	void part1() {
		ColorImage image = loadColorPPM("peppers.ppm");
		GrayscaleImage grayscale = convertToGrayscale(image);
		saveGrayscalePPM("grayscale_peppers.ppm", grayscale);
	}
}


int main() {
	lab5::part1();	
}