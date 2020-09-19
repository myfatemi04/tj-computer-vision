#include<iostream>
#include<fstream>
#include<cstdlib>
#include<math.h>
#include<time.h>

using namespace std;

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

        // Equivalent to (dy / dx) - 1
        int error = dy - dx;

        for (int i = x1; i <= x2; i++) {
            pixels[i][j] = color;
            
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
            
            if (error > 0 && x2 > x1) {
                j += 1;

                error -= dy;
            } else if (error < 0 && x2 < x1) {
                j -= 1;

                error += dy;
            }

            error += dx;
        }
    }
}

void drawCircle(int*** pixels, int center_x, int center_y, int radius, int* color) {
    cout << "Drawing a circle centered at (" << center_x << ", " << center_y << ") with radius " << radius;

    int x, y, max_x, y2, y2_new, two_y;

    // we draw the pixels from the top until 45 degrees clockwise from the top
    max_x = (int) (radius * 0.70710678);

    // starts with the topmost point
    x = 0;
    y = radius;

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

        // Set all symmetric points to 1
        pixels[x + center_x][y + center_y] = color;
        pixels[-x + center_x][y + center_y] = color;
        pixels[x + center_x][-y + center_y] = color;
        pixels[-x + center_x][-y + center_y] = color;

        pixels[y + center_x][x + center_y] = color;
        pixels[-y + center_x][x + center_y] = color;
        pixels[y + center_x][-x + center_y] = color;
        pixels[-y + center_x][-x + center_y] = color;

        y2_new -= 2 * x - 3;

        x += 1;
    }
}

double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

double* findIncircle(double ax, double ay, double bx, double by, double cx, double cy) {
    double A = distance(bx, by, cx, cy);
    double B = distance(ax, ay, cx, cy);
    double C = distance(bx, by, ax, ay);

    double p = A + B + C;
    double s = p/2;
    double area = sqrt((s - A) * (s - B) * (s - C) / s);
    double inradius = area;

    double ox = (A * ax + B * bx + C * cx) / p;
    double oy = (A * ay + B * by + C * cy) / p;

    double* incircle = new double[3];

    incircle[0] = ox;
    incircle[1] = oy;
    incircle[2] = inradius;

    cout << ox << " " << oy << " " << inradius << "\n";

    return incircle;
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

double getRandom() {
    return rand() / ((double) RAND_MAX);
}

double scale(double value, int max) {
    return (int) (value * max);
}

int main() {
    // set up random seed
    srand(time(nullptr));

    int width = 200;
    int height = 200;

    // initialize the image array
    int*** myArray = new int**[width];
    initializeArray(myArray, width, height);

    // These should be pointers
    int* WHITE = new int[3];
    WHITE[0] = WHITE[1] = WHITE[2] = 1;

    int* RED = new int[3];
    RED[0] = 1;
    RED[1] = RED[2] = 0;

    int* GREEN = new int[3];
    GREEN[0] = GREEN[2] = 0;
    GREEN[1] = 1;

    int* BLUE = new int[3];
    BLUE[0] = BLUE[1] = 0;
    BLUE[2] = 1;

    int* PURPLE = new int[3];
    PURPLE[0] = PURPLE[2] = 1;
    PURPLE[1] = 0;

    // Initialize random triangle points
    // Internally, these will be points inside a unit square
    double x1, y1, x2, y2, x3, y3;
    int x1_scale, y1_scale, x2_scale, y2_scale, x3_scale, y3_scale;
    x1 = getRandom();
    y1 = getRandom();
    x2 = getRandom();
    y2 = getRandom();
    x3 = getRandom();
    y3 = getRandom();

    x1_scale = scale(x1, width);
    y1_scale = scale(y1, height);
    x2_scale = scale(x2, width);
    y2_scale = scale(y2, height);
    x3_scale = scale(x3, width);
    y3_scale = scale(y3, height);

    // Draw the triangle
    drawBresenhams(myArray, x1_scale, y1_scale, x2_scale, y2_scale, RED);
    drawBresenhams(myArray, x1_scale, y1_scale, x3_scale, y3_scale, GREEN);
    drawBresenhams(myArray, x2_scale, y2_scale, x3_scale, y3_scale, BLUE);

    myArray[x1_scale][y1_scale] = PURPLE;
    myArray[x2_scale][y2_scale] = PURPLE;
    myArray[x3_scale][y3_scale] = PURPLE;

    // Draw incircle
    double* incircle = findIncircle(x1, y1, x2, y2, x3, y3);
    drawCircle(myArray,
        scale(incircle[0], width), // x
        scale(incircle[1], height), // y
        scale(incircle[2], width), // radius
        WHITE);

    // save the image
    ofstream fileout("output.ppm");
    fileout << "P3 " << width << " " << height << " 1" << "\n";

    printArray(myArray, fileout, width, height);

    fileout.close();
}