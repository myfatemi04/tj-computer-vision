#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

int main(int argc, char** argv) {
	std::cout << "Hello!\n";

	if (argc != 2) {
		std::cout << "usage: DisplayImage <Image_Path>\n";
		return -1;
	}

	std::cout << "Reading image from " << argv[1] << '\n';

	Mat image;
	image = imread(argv[1], ImreadModes::IMREAD_COLOR);

	if (!image.data) {
		printf("No image data \n");
		return -1;
	}
	/*
	//next 2 lines would be to create a window and display in the window the image
	//you can uncomment them if you install openCV on your computer
	namedWindow("Display Image", WINDOW_AUTOSIZE );
	imshow("Display Image", image);
	*/
	cout << "Let's create a copy now..." << endl;
	//now lets save the Mat in a new file
	imwrite("./copy.jpg", image);
	
	//now lets create the grayscale image
	Mat gray;
	
	cout << "Let's create a grayscale now..." << endl;
	gray = imread( argv[1], IMREAD_GRAYSCALE );
	imwrite("./gray.jpg", gray);	

	waitKey(0);

	return 0;
}