// l071.cpp
// Michael Fatemi

#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

int main(int argc, char** argv) {
	if (argc != 2) {
		std::cout << "usage: l071 <path>\n";
		std::cout << "path: the source image for detecting coins\n";
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

	const int GAUSSIAN_KERNEL_SIZE = 5;
	
	cv::GaussianBlur(
		gray, // input
		gray, // output
		cv::Size(GAUSSIAN_KERNEL_SIZE, GAUSSIAN_KERNEL_SIZE), // kernel size
		3, // sigma X
		3, // sigma Y
		cv::BORDER_DEFAULT // border style
	);

	const int MIN_RADIUS = 75;
	const int MAX_RADIUS = 150;
	const int MIN_DISTANCE = MIN_RADIUS;

	const int MAX_POSSIBLE_EDGE_GRADIENT = 1016;

	const int EDGE_DETECTION_THRESHOLD = MAX_POSSIBLE_EDGE_GRADIENT * 0.77 / GAUSSIAN_KERNEL_SIZE;
	const int ACCUMULATOR_THRESHOLD = 40;

	const double ACCUMULATOR_BLOCKINESS = 1.5;

  std::vector<cv::Vec3f> circles;

	// cv: Detect circles in [gray] at least [MIN_DISTANCE] apart with a radius from [MIN_RADIUS] to [MAX_RADIUS] and output them to [circles]
	
	cv::Mat edges;
	// cv: Detect edges in (gray) with thresholds (EDGE_DETECTION_THRESHOLD/2, EDGE_DETECTION_THRESHOLD)
	// Replicate the same style of edge detection used by Hough Circles
	cv::Canny(gray, edges, EDGE_DETECTION_THRESHOLD / 2, EDGE_DETECTION_THRESHOLD);

	cv::imwrite("edges.jpg", edges);

	cv::HoughCircles(
		gray, // image
		circles, // output
		cv::HOUGH_GRADIENT, // method
		ACCUMULATOR_BLOCKINESS, // image resolution / accumulator resolution
		MIN_DISTANCE,  // minimum distance
		EDGE_DETECTION_THRESHOLD, // edge detection threshold
		ACCUMULATOR_THRESHOLD,  // accumulator threshold
		MIN_RADIUS, // minimum radius
		MAX_RADIUS  // maximum radius
	);

	for (const auto& circle : circles) {
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
