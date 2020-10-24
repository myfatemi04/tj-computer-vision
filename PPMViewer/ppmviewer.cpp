#include <windows.h>
// #include <glad/glad.h>
#include <gl/GL.h>
#include <gl/GLU.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <fstream>

using namespace std;

class InvalidPixelException: public exception {
    private:
        int x, y, value;
        const char* filetype;
    public:
        InvalidPixelException(int x_, int y_, int value_, const char* filetype_) {
            x = x_;
            y = y_;
            value = value_;
            filetype = filetype_;
        }
        int getX() {
            return x;
        }
        int getY() {
            return y;
        }
        int getValue() {
            return value;
        }
        const char* getFiletype() {
            return filetype;
        }
};

class PPMImage {
    private:
        unsigned char*** pixels;
        int width, height, channels;
        int type;

    public:
        PPMImage(unsigned char*** pixels, int width, int height, int channels, int type) {
            this -> pixels = pixels;
            this -> width = width;
            this -> height = height;
            this -> channels = channels;
            this -> type = type;
        }

        unsigned char*** getPixels() {
            return pixels;
        }

        int getWidth() const {
            return width;
        }

        int getHeight() const {
            return height;
        }
};

bool DEBUG = true;

bool validByte(int x) {
    return (x >= 0) && (x <= 255);
}

PPMImage readPPM(istream& in) {
    in.ignore(1); // "P"
    int type, width, height, grayLevel = 1, channels = 1;
    in >> type;
    in >> width;
    in >> height;

    if (DEBUG) {
        cout << "Width and height: " << width << " x " << height << "\n";
    }

    // This is a PBM (binary image)
    if (type == 1) {
        if (DEBUG) {
            cout << "Type: BINARY\n";
        }
    } else if (type == 2) {
        if (DEBUG) {
            cout << "Type: GRAYSCALE\n";
        }
        in >> grayLevel;
    } else if (type == 3) {
        if (DEBUG) {
            cout << "Type: COLOR\n";
        }
        in >> grayLevel;
        channels = 3;
    }

    unsigned char*** pixels = new unsigned char**[height];
    for (int y = 0; y < height; y++) {
        unsigned char** row = new unsigned char*[width];
        for (int x = 0; x < width; x++) {
            // Every pixel will be an RGB pixel, even if the PPM is binary
            unsigned char* pixel = new unsigned char[3];

            if (channels == 1) {
                // Read one integer
                int v = 0;
                in >> v;

                if (type == 1) {
                    if (v != 0 || v != 1) {
                        throw InvalidPixelException(x, y, v, "binary");
                    }

                    // v = 0 or 1, so we scale to 0-255.
                    v = v * 255;
                } else if (type == 2) {
                    if (!validByte(v)) {
                        throw InvalidPixelException(x, y, v, "grayscale");
                    }
                }

                pixel[0] = pixel[1] = pixel[2] = v;
            } else {
                int r, g, b;
                in >> r >> g >> b;

                if (!validByte(r)) {
                    throw InvalidPixelException(x, y, r, "color");
                }

                if (!validByte(g)) {
                    throw InvalidPixelException(x, y, g, "color");
                }

                if (!validByte(b)) {
                    throw InvalidPixelException(x, y, b, "color");
                }

                pixel[0] = r;
                pixel[1] = g;
                pixel[2] = b;
            }
            row[x] = pixel;
        }
        pixels[y] = row;
    }

    return PPMImage(pixels, width, height, channels, type);
}

void setSize(int width, int height) {
    // (lower left), (upper right)
    glViewport(0, height, width, 0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluOrtho2D(0, width, height, 0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

int main(int argc, char *argv[]) {
    if (argc == 1) {
        cout << "Usage: ppmviewer <filename>" << endl;
        return 0;
    }

    char* filename = argv[1];
    
    if (!glfwInit()) {
        cerr << "Failed to initialize GLFW\n";
        return 1;
    }

    /* If you include these it will fail
    // I'm using Intel GPU but I don't care enough to use glfwWindowHint

    // idk what this does but it was in all the tutorials
    glfwWindowHint(GLFW_SAMPLES, 4);
    // Hint that we are using GLFW 3.1
    glfwWindowHint(GLFW_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_VERSION_MINOR, 0);
    // Use OpenGL core profile
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // To make MacOS work
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    */

    GLFWwindow* window;
    window = glfwCreateWindow(800, 800, "PPM", NULL, NULL);

    if (window == NULL) {
        cerr << "Failed to open GLFW window.\n";
        return 1;
    }
    
    glfwMakeContextCurrent(window);

    // ifstream imageFile(filename);
    cout << "Reading from " << filename << "\n";

    ifstream imageFile(filename);

    if (imageFile.fail()) {
        cerr << "Image file failed to load: " << filename << endl;
        return 1;
    }

    PPMImage myImage = readPPM(imageFile);

    glClear(GL_COLOR_BUFFER_BIT); // Clear the display
    glRasterPos2i(0, 0); // Start drawing the image from (0, 0)
    glDrawPixels(myImage.getWidth(), myImage.getHeight(), GL_RGB, GL_UNSIGNED_BYTE, myImage.getPixels());
    
    while(!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);

        glfwSwapBuffers(window);
        glfwPollEvents();    
    }
}