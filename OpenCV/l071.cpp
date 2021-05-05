#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

int main(int argc, char** argv) {
	if (argc != 2) {
		std::cout << "usage: l071 <path>\n";
		std::cout << "path: the path of the image to run Hough Circles on\n";
		return -1;
	}

	std::cout << "Reading image from " << argv[1] << '\n';

	cv::Mat image = cv::imread(argv[1], cv::ImreadModes::IMREAD_COLOR);

	if (image.empty()) {
		std::cout << "No image data\n";
		return -1;
	}

	cv::Mat gray;
	cv::cvtColor(image, gray, cv::COLOR_BGR2GRAY);
	
	cv::medianBlur(gray, gray, 5);

	const int MIN_RADIUS = 50;
	const int MAX_RADIUS = 300;
	const int MIN_DISTANCE = MIN_RADIUS;

	const int EDGE_DETECTION_THRESHOLD = 100;
	const int ACCUMULATOR_THRESHOLD = 30;

	const int DP = 1 / 1;

  std::vector<cv::Vec3f> circles;

	// cv: Detect circles in [gray] at least [MIN_DISTANCE] apart with a radius from [MIN_RADIUS] to [MAX_RADIUS] and output them to [circles]

	cv::HoughCircles(
		gray, // image
		circles, // output
		cv::HOUGH_GRADIENT, // method
		DP, // dp
		MIN_DISTANCE,  // minimum distance
		EDGE_DETECTION_THRESHOLD, // edge detection threshold
		ACCUMULATOR_THRESHOLD,  // accumulator threshold
		MIN_RADIUS, // minimum radius
		MAX_RADIUS  // maximum radius
	);

	for(const auto& circle : circles) {
		int x = circle[0];
		int y = circle[1];
		int radius = circle[2];

		auto center = cv::Point(x, y);
		cv::circle(image, center, 1, cv::Scalar(0, 100, 100), 3, cv::LINE_AA);
		cv::circle(image, center, radius, cv::Scalar(255, 0, 255), 3, cv::LINE_AA);
	}

	cv::imwrite("imagec.jpg", image);

	cv::waitKey(0);

	return 0;
}
