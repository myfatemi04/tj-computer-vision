#include <iostream>
#include <fstream>
using namespace std;

void save_image(int ***pixels, int width, int height, int max, const char* filename) {
   ofstream file;
   file.open(filename);
   file << "P3 " << width << " " << height << " " << max << "\n";
   
   for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
         for (int c = 0; c < 3; c++) {
            int value = pixels[y][x][c];
            file << value << " ";
         }
      }
      file << "\n";
   }
   
   file.close();
}

int main() {
   int ***pixels = new int[1][1][3];
   pixels[0][0][1] = 1;
   save_image(pixels, 1, 1, 1, "test.ppm");
   
   return 0;
}