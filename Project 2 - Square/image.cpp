#include<iostream>
#include<vector>
#include "geometry.cpp"

void swap(int* a, int* b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

class Image {
    private:
        int const*** pixels;
        int width, height, channels;

    public:
        int getWidth() {
            return width;
        }

        int getHeight() {
            return height;
        }

        int getChannels() {
            return channels;
        }

        Image(int width_, int height_, int channels_, const int* bgColor):  
            width(width_),
            height(height_),
            channels(channels_) {
            pixels = new int const**[height];

            // allocate arrays for all pixels
            for (int y = 0; y < height; y++) {
                pixels[y] = new int const*[width];
                for (int x = 0; x < width; x++) {
                    pixels[y][x] = bgColor;
                }
            }

        }

        void putPixel(int x, int y, const int* color) {
            if (x < 0 || x >= width) return;
            if (y < 0 || y >= height) return;

            pixels[y][x] = color;
        }

        const int* getPixel(int x, int y) {
            return pixels[y][x];
        }

        void drawLineSegment(Point a, Point b, const int* color) {
            int x1 = (int)a.getX();
            int y1 = (int)a.getY();
            int x2 = (int)b.getX();
            int y2 = (int)b.getY();

            // printf("Drawing from (%d, %d) to (%d, %d)\n", x1, y1, x2, y2);

            int dx = x2 - x1;
            int dy = y2 - y1;

            bool isXMajor = abs(dx) > abs(dy);

            if (isXMajor) {
                if (x1 > x2) {
                    swap(&x1, &x2);
                    swap(&y1, &y2);
                    dx = -dx; // we swapped points, so these must be updated
                    dy = -dy;
                }

                int j = y1;
                int error = dy; // Equivalent to (dy / dx) - 1

                if (dy > 0) {
                    error -= dx;
                } else {
                    error += dx;
                }

                for (int i = x1; i <= x2; i++) {
                    putPixel(i, j, color);
                    
                    // Error is greater than 0, so we know the line
                    // has moved up one pixel
                    if (error > 0 && y2 > y1) {
                        j += 1;

                        error -= dx;
                    } else if (error < 0 && y2 < y1) {
                        j -= 1;

                        error += dx;
                    }

                    error += dy;
                }

                putPixel(x2, j, color);
            } else {
                if (y1 > y2) {
                    swap(&x1, &x2);
                    swap(&y1, &y2);
                    dx = -dx; // we swapped points, so these must be updated
                    dy = -dy;
                }

                int j = x1;
                int error = dx; // Equivalent to (dx / dy) - 1

                if (dx > 0) {
                    error -= dy;
                } else {
                    error += dy;
                }

                for (int i = y1; i <= y2; i++) {
                    putPixel(j, i, color);
                    
                    if (error > 0 && x2 > x1) {
                        j += 1;

                        error -= dy;
                    } else if (error < 0 && x2 < x1) {
                        j -= 1;

                        error += dy;
                    }

                    error += dx;
                }

                putPixel(j, y2, color);
            }
        }

        void drawLineSegment(LineSegment segment, const int* color) {
            drawLineSegment(segment.getA(), segment.getB(), color);
        }

        void drawCircle(Circle circle, const int* color) {
            double radius = circle.getRadius();
            int center_x = (int) circle.getCenter().getX();
            int center_y = (int) circle.getCenter().getY();

            int x, y, y2, y2_new, two_y;

            // starts with the topmost point
            x = 0;
            y = (int)(radius + 0.5); // round up

            y2 = y * y;
            y2_new = y2;

            two_y = 2 * y - 1;

            while (y >= x) {
                // when X increases, see if the Y value should decrease.
                if ((y2 - y2_new) >= two_y) {
                    y2 -= two_y;

                    // decrease Y and 2Y
                    y -= 1;
                    two_y -= 2;
                }

                putPixel(x + center_x, y + center_y, color);
                putPixel(x + center_x, -y + center_y, color);
                putPixel(-x + center_x, y + center_y, color);
                putPixel(-x + center_x, -y + center_y, color);

                putPixel(y + center_x, x + center_y, color);
                putPixel(y + center_x, -x + center_y, color);
                putPixel(-y + center_x, x + center_y, color);
                putPixel(-y + center_x, -x + center_y, color);

                y2_new -= 2 * x - 3;

                x += 1;
            }
        }

        void drawPolygon(Polygon polygon, const int* color, bool drawCircles = false) {
            std::vector<LineSegment> sides = polygon.getSides();

            for (int i = 0; i < polygon.length(); i++) {
                drawLineSegment(polygon.getSide(i), color);
            }

            if (drawCircles) {
                for (int i = 0; i < polygon.length(); i++) {
                    drawCircle(Circle(polygon.getPoint(i), 2), color);
                }
            }
        }

        void saveImage(std::ostream& out) {
            out << "P3 " << width << " " << height << " 1\n";

            int x, y, c;

            for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                    for (c = 0; c < channels; c++) {
                        out << pixels[y][x][c] << " ";
                    }
                }

                out << "\n";
            }
        }
};