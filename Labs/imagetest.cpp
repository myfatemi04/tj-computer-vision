#include "image.cpp"
#include<iostream>
#include<fstream>

using namespace std;

int main() {
    int * WHITE = new int[3];
    WHITE[0] = WHITE[1] = WHITE[2] = 1;

    Image x(200, 200, 3, WHITE);
    
    x.drawLine(10, 10, 50, 50, WHITE);

    ofstream fileout("test.ppm");

    x.saveImage(fileout);

    fileout.close();

    cout << "Saved image to test.ppm" << endl;
}

