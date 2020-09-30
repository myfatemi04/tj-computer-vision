#include<iostream>
#include<fstream>
#include<cstdlib>
#include<math.h>
#include<time.h>

using namespace std;

class Circle {
    public:
        double x, y, radius;
};

void putPixel(int*** pixels, int* color, int x, int y, int width, int height) {
    if (x < 0 || x >= width) return;
    if (y < 0 || y >= height) return;

    pixels[x][y] = color;
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

        int j = y1;

        // Equivalent to (dy / dx) - 1
        int error = dy - dx;

        for (int i = x1; i <= x2; i++) {
            putPixel(pixels, color, i, j, 800, 800);
            
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

        int j = x1;
        int error = dx - dy;

        for (int i = y1; i <= y2; i++) {
            putPixel(pixels, color, j, i, 800, 800);
            
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

void drawCircle(int*** pixels, int center_x, int center_y, double radius, int* color) {
    cout << "Drawing a circle centered at (" << center_x << ", " << center_y << ") with radius " << radius << "\n";

    int x, y, y2, y2_new, two_y;

    // starts with the topmost point
    x = 0;
    y = round(radius);

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
        int width = 800, height = 800;

        putPixel(pixels, color, x + center_x, y + center_y, width, height);
        putPixel(pixels, color, x + center_x, -y + center_y, width, height);
        putPixel(pixels, color, -x + center_x, y + center_y, width, height);
        putPixel(pixels, color, -x + center_x, -y + center_y, width, height);

        putPixel(pixels, color, y + center_x, x + center_y, width, height);
        putPixel(pixels, color, y + center_x, -x + center_y, width, height);
        putPixel(pixels, color, -y + center_x, x + center_y, width, height);
        putPixel(pixels, color, -y + center_x, -x + center_y, width, height);

        y2_new -= 2 * x - 3;

        x += 1;
    }
}

double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

Circle findIncircle(double ax, double ay, double bx, double by, double cx, double cy) {
    double A = distance(bx, by, cx, cy);
    double B = distance(ax, ay, cx, cy);
    double C = distance(ax, ay, bx, by);

    double p = A + B + C;
    double s = p / 2;
    double area = sqrt(((s - A) * (s - B) * (s - C)) / s);
    double inradius = area;

    Circle myCircle;

    double ox = (A * ax + B * bx + C * cx) / p;
    double oy = (A * ay + B * by + C * cy) / p;

    myCircle.x = ox;
    myCircle.y = oy;
    myCircle.radius = inradius;

    return myCircle;
}

double* solveLinear(double a1, double b1, double c1, double a2, double b2, double c2) {
    // matrix is:
    /*
    | a1 b1 | = | c1 |
    | a2 b2 |   | c2 |

    We will use Cramer's rule to find the (x, y) of the intersection.
    */

    double det = a1 * b2 - b1 * a2;

    if (det == 0) {
        throw 0;
    }

    double detX = c1 * b2 - b1 * c2;
    double detY = a1 * c2 - c1 * a2;

    double* solution = new double[2] { detX / det, detY / det };

    return solution;
}

// This function uses linear equations to find the circumcircle
// Circumcircle = intersection of perpendicular bisectors
Circle findCircumcircle(double ax, double ay, double bx, double by, double cx, double cy) {
    // find midpoints
    double ab_x = (ax + bx) / 2;
    double ab_y = (ay + by) / 2;
    double ac_x = (ax + cx) / 2;
    double ac_y = (ay + cy) / 2;

    /*
    equations are ax + by = c
    how to find a, b, and c?
    y = mx + b
    m = slope_numer/slope_denom
    y = slope_numer/slope_denom * x + b
    slope_denom(y) = slope_numer(x) + b
    -slope_numer(x) + slope_denom(y) = b
    now, substitute a point for (x, y) to find b. we can just use (x1, y1)
    b = -slope_numer(x1) + slope_denom(y1)
    now, we know that in ax + by = c:
    a = -slope_numer
    b = slope_denom
    c = -a(x1) + b(y1)
    */
    
    // cout << "(" << ax << ", " << ay << ") (" << bx << ", " << by << ") (" << cx << ", " << cy << ")\n";

    // find a, b, and c of first side
    double ab_slope_numer = -(bx - ax);
    double ab_slope_denom = (by - ay);

    // cout << "AB slope numer: " << ab_slope_numer << endl;
    // cout << "AB slope denom: " << ab_slope_denom << endl;

    double ab_a = -ab_slope_numer;
    double ab_b = ab_slope_denom;
    double ab_c = ab_a * ab_x + ab_b * ab_y;

    // find a, b, and c of second side
    double ac_slope_numer = -(cx - ax);
    double ac_slope_denom = (cy - ay);

    double ac_a = -ac_slope_numer;
    double ac_b = ac_slope_denom;
    double ac_c = ac_a * ac_x + ac_b * ac_y;

    // solve the linear equation for the center
    double circumcenter_x = 0;
    double circumcenter_y = 0;

    try {
        double* circumcenter = solveLinear(ab_a, ab_b, ab_c, ac_a, ac_b, ac_c);
        
        circumcenter_x = circumcenter[0];
        circumcenter_y = circumcenter[1];
    } catch (int i) {
        if (i == 0) {
            // this means there is no solution (the lines are parallel)
            // the only way for the lines to be parallel in a triangle is
            // if the lines are on top of each other
            
            // so, we can just use a midpoint
            circumcenter_x = ab_x;
            circumcenter_y = ab_y;
        }
    }

    // Now, we must find the radius. We will simply find the distance
    // from the circumcenter to one of the points.
    double radius = distance(circumcenter_x, circumcenter_y, ax, ay);

    Circle myCircle{circumcenter_x, circumcenter_y, radius};

    return myCircle;
}

Circle find9PointCircle(double ax, double ay, double bx, double by, double cx, double cy) {
    // The 9 point circle goes through the midpoints of each side
    // So, it is the circumcircle of the midpoints of each side

    // find midpoints
    double ab_x = (ax + bx) / 2;
    double ab_y = (ay + by) / 2;
    double ac_x = (ax + cx) / 2;
    double ac_y = (ay + cy) / 2;
    double bc_x = (bx + cx) / 2;
    double bc_y = (by + cy) / 2;

    return findCircumcircle(ab_x, ab_y, ac_x, ac_y, bc_x, bc_y);
}

double* findCentroid(double ax, double ay, double bx, double by, double cx, double cy) {
    double* point = new double[2];
    point[0] = (ax + bx + cx) / 3;
    point[1] = (ay + by + cy) / 3;

    return point;
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

    int width = 800;
    int height = 800;

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
    Circle incircle = findIncircle(x1, y1, x2, y2, x3, y3);
    drawCircle(myArray,
        scale(incircle.x, width), // x
        scale(incircle.y, height), // y
        incircle.radius * width, // radius
        WHITE);

    // Draw circumcircle
    Circle circumcircle = findCircumcircle(x1, y1, x2, y2, x3, y3);
    drawCircle(myArray,
        scale(circumcircle.x, width), // x
        scale(circumcircle.y, height), // y
        circumcircle.radius * width, // radius
        WHITE);

    // Draw 9 point circle
    Circle ninePointCircle = find9PointCircle(x1, y1, x2, y2, x3, y3);
    drawCircle(myArray,
        scale(ninePointCircle.x, width), // x
        scale(ninePointCircle.y, height), // y
        ninePointCircle.radius * width, // radius
        WHITE);

    // Find centroid
    double* centroid = findCentroid(x1, y1, x2, y2, x3, y3);

    // The euler line is the line between the centroid and circumcenter
    // When we draw it, we must extend it to the edges of the screen
    double lineSegmentLength = distance(centroid[0], centroid[1], circumcircle.x, circumcircle.y);

    // we want to scale this line segment to be length sqrt2/2
    double increase = sqrt(2) / lineSegmentLength;

    double elAX = centroid[0];
    double elAY = centroid[1];
    double elBX = circumcircle.x;
    double elBY = circumcircle.y;

    double nelAX = elAX - (elBX - elAX) * increase;
    double nelAY = elAY - (elBY - elAY) * increase;
    double nelBX = elBX + (elBX - elAX) * increase;
    double nelBY = elBY + (elBY - elAY) * increase;
    
    drawBresenhams(myArray,
        scale(nelAX, width),
        scale(nelAY, height),
        scale(nelBX, width),
        scale(nelBY, height),
        PURPLE);

    // save the image
    ofstream fileout("triangle.ppm");
    fileout << "P3 " << width << " " << height << " 1" << "\n";

    printArray(myArray, fileout, width, height);

    fileout.close();
}