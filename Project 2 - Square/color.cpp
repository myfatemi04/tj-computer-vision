#pragma once

class Color {
    private:
        int maxValues;
        int* values;
        int channelCount;

    public:
        Color(int* values, int channelCount = 3, int maxValues = 1) {
            this -> values = values;
            this -> channelCount = channelCount;
            this -> maxValues = maxValues;
        }

        int* getPixel() {
            int* pixel = new int[channelCount];
            for (int i = 0; i < channelCount; i++) {
                pixel[i] = values[i];
            }
            return pixel;
        }
};

namespace Colors {
    const int* WHITE = new int[3] { 1, 1, 1 };
    const int* BLACK = new int[3] { 0, 0, 0 };
    const int* RED = new int[3] { 1, 0, 0 };
    const int* GREEN = new int[3] { 0, 1, 0 };
    const int* BLUE = new int[3] { 0, 0, 1 };
};