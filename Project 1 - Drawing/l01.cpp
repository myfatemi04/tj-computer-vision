#include<iostream>
#include<fstream>
#include<cstdlib>
#include<math.h>
#include<time.h>

using namespace std;

void swap(int* a, int* b) {
    int buf = *a;
    a = b;
    b = a;
}

void drawBresenhams(int*** pixels, int x1, int y1, int x2, int y2, int* color) {
    cout << "Drawing line connecting (" << x1 << ", " << y1 << ") and (" << x2 << ", " << y2 << ").\n";

    int dx = x2 - x1;
    int dy = y2 - y1;

    bool isXMajor = abs(dx) > abs(dy);

    if (isXMajor) {
        if (x1 > x2) {
            swap(x1, x2);
            swap(y1, y2);
            dx = -dx; // we swapped points, so these must be updated
            dy = -dy;
        }

        cout << "X Major\n";

        int j = y1;
        int error = dy - dx;

        for (int i = x1; i <= x2; i++) {
            pixels[i][j] = color;
            
            if (error >= 0) {
                if (y2 > y1) {
                    j += 1;
                } else if (y2 < y1) {
                    j -= 1;
                }

                error -= dx;
            }

            error += dy;
        }
    } else {
        if (y1 > y2) {
            swap(x1, x2);
            swap(y1, y2);
            dx = -dx; // we swapped points, so these must be updated
            dy = -dy;
        }

        cout << "Y Major\n";

        int j = x1;
        int error = dx - dy;

        for (int i = y1; i <= y2; i++) {
            pixels[j][i] = color;
            
            if (error >= 0) {
                if (x2 > x1) {
                    j += 1;
                } else if (x2 < x1) {
                    j -= 1;
                }

                error -= dy;
            }

            error += dx;
        }
    }
}

void drawCircle(int*** pixels, int center_x, int center_y, int radius, int* color) {
    int x, y, max_x, y2, y2_new, two_y;

    // we draw the pixels from the top until 45 degrees clockwise from the top
    max_x = (int) (radius * 0.70710678);

    // starts with the topmost point
    x = 0;
    y = radius;

    y2 = y * y;
    y2_new = y2;

    two_y = 2 * y - 1;

    for (x = 0; x < max_x; x++) {
        // when X increases, see if the Y value should decrease.
        if ((y2 - y2_new) >= two_y) {
            y2 -= two_y;

            // decrease Y and 2Y
            y -= 1;
            two_y -= 2;
        }

        // Set all symmetric points to 1
        pixels[x][y] = color;
        pixels[-x][y] = color;
        pixels[x][-y] = color;
        pixels[-x][-y] = color;

        pixels[y][x] = color;
        pixels[-y][x] = color;
        pixels[y][-x] = color;
        pixels[-y][-x] = color;

        y2_new -= 2 * x - 3;
    }
}

void printArray(int*** array, ostream& out, int width, int height, int channels = 3) {
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int z = 0; z < channels; z++) {
                out << array[x][y][z] << " ";
            }
        }

        out << "\n";
    }
}

void initializeArray(int*** array, int width, int height, int channels = 3) {
    for (int i = 0; i < width; i++) {
        array[i] = new int*[height];
        for (int j = 0; j < height; j++) {
            array[i][j] = new int[channels];
            for (int k = 0; k < channels; k++) {
                array[i][j][k] = 0;
            }
        }
    }
}

int getRandom(int min, int max) {
    int range = max - min;
    return min + range * (rand() / ((double) RAND_MAX));
}

int main() {
    // set up random seed
    srand(time(nullptr));

    int width = 800;
    int height = 800;

    int*** myArray = new int**[width];

    int* WHITE = new int[3];
    WHITE[0] = 1;
    WHITE[1] = 1;
    WHITE[2] = 1;

    ofstream fileout("output.ppm");
    fileout << "P3 " << width << " " << height << " 1" << "\n";

    initializeArray(myArray, width, height);

    int x1, y1, x2, y2;
    x1 = getRandom(0, width - 1);
    y1 = getRandom(0, height - 1);
    x2 = getRandom(0, width - 1);
    y2 = getRandom(0, height - 1);

    drawBresenhams(myArray, x1, y1, x2, y2, WHITE);

    printArray(myArray, fileout, width, height);

    fileout.close();
}