// l071.cpp
// Michael Fatemi

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

typedef std::vector<cv::Vec3f> Circles;

double _safeAbs(double a) {
	return a * (2 * (a > 0) - 1);
}

double absoluteError(double a, double b) {
		return _safeAbs(a - b);
	}

namespace coins {
	enum CoinType {
		COIN_TYPE_PENNY = 0,
		COIN_TYPE_NICKEL = 1,
		COIN_TYPE_DIME = 2,
		COIN_TYPE_QUARTER = 3,
		COIN_TYPE_HALF_DOLLAR = 4,
		COIN_TYPE_SILVER_DOLLAR = 5
	};

	const double PENNY_RADIUS_MM = 19.05;

	const double* COIN_SIZE_RELATIVE_TO_PENNY = new double[6] {
		1,                       // penny
		21.21 / PENNY_RADIUS_MM, // nickel
		17.91 / PENNY_RADIUS_MM, // dime
		24.26 / PENNY_RADIUS_MM, // quarter
		30.61 / PENNY_RADIUS_MM, // half dollar
		38.1 / PENNY_RADIUS_MM   // silver dollar
	};

	const CoinType* COIN_TYPES = new CoinType[6] {
		COIN_TYPE_PENNY,
		COIN_TYPE_NICKEL,
		COIN_TYPE_DIME,
		COIN_TYPE_QUARTER,
		COIN_TYPE_HALF_DOLLAR,
		COIN_TYPE_SILVER_DOLLAR
	};

	const int COIN_COLORS[6][3] = {
		{ 255,   0,   0 }, // penny: red
		{ 255,   0, 255 }, // nickel: purple
		{   0,   0, 255 }, // dime: blue
		{   0, 255,   0 }, // quarter: green
		{   0,   0,   0 }, // half dollar: unknown
		{ 255, 255,   0 }  // silver dollar: yellow
	};

	const int COIN_VALUES[6] = {
		1, 5, 10, 25, 50, 100
	};

	CoinType estimateCoinType(double unknownCoinRadius, double pennyRadius) {
		double radiusRelativeToPenny = unknownCoinRadius / pennyRadius;

		CoinType leastErrorCoinType = COIN_TYPE_PENNY;
		double leastErrorCoinTypeError = absoluteError(1, unknownCoinRadius);
		
		for (int i = 0; i < 6; i++) {
			double error = absoluteError(radiusRelativeToPenny, COIN_SIZE_RELATIVE_TO_PENNY[i]);
			if (error < leastErrorCoinTypeError) {
				leastErrorCoinTypeError = error;
				leastErrorCoinType = COIN_TYPES[i];
			}
		}

		return leastErrorCoinType;
	}

	// Estimates by assuming that the smallest radius is that of a dime
	double findPennyRadiusDimeRadiusFallback(const Circles& circles) {
		std::cout << "Finding penny radius went to fallback\n";

		int dimeRadius = -1;
		// circle is [x, y, radius]
		for (const auto& circle : circles) {
			int radius = circle[2];
			if (dimeRadius == -1 || radius < dimeRadius) {
				dimeRadius = radius;
			}
		}

		if (dimeRadius == -1) {
			return -1;
		}

		return dimeRadius / COIN_SIZE_RELATIVE_TO_PENNY[COIN_TYPE_DIME];
	}

	double findPennyRadius(const Circles& circles, const cv::Mat& reference) {
		std::cout << "Finding penny radius\n";

		double pennyRadiiCumulativeTotal = 0;
		int pennyCount = 0;
		for (const auto& circle : circles) {
			int x = circle[0];
			int y = circle[1];
			int radius = circle[2];
			auto center = cv::Point(x, y);

			// cv: Get points of ellipse centered at (center) with axes (radius/2, radius/2) [and resolution 1]
			std::vector<cv::Point> points;
			cv::ellipse2Poly(center, cv::Size(radius / 2, radius / 2), 0, 0, 360, 1, points);

			int redPixelCount = 0;
			int totalPixelCount = points.size();
			// for each (location) of (points)
			for (const auto& point : points) {
				auto color = reference.at<cv::Vec3b>(point);
				// [blue, green, red] = color
				int blue = color[0];
				int green = color[1];
				int red = color[2];
				double averageColor = (blue + green + red) / 3.0;
				if (red > (averageColor * 1.3)) {
					redPixelCount++;
				}
			}
			if (redPixelCount > totalPixelCount * 0.3) {
				pennyCount++;
				pennyRadiiCumulativeTotal += radius;
			}
		}
		
		if (pennyCount == 0) {
			// Assume the smallest coin is a dime, and extrapolate the penny
			// radius from there.
			return findPennyRadiusDimeRadiusFallback(circles);
		} else {
			// Return the average radius of the found pennies
			return pennyRadiiCumulativeTotal / pennyCount;
		}
	}
}

std::string formatMoney(int cents) {
	int dollarCount = cents / 100;
	int centCount = cents % 100;

	std::stringstream ss;
	ss << '$' << dollarCount << '.' << (centCount < 10 ? "0" : "") << centCount;
	
	return ss.str();
}

void writeResults(std::ostream& out, const int* quantities, int cents) {
	out << "Result:\n";
	out << quantities[coins::COIN_TYPE_PENNY] << " pennies\n";
	out << quantities[coins::COIN_TYPE_NICKEL] << " nickels\n";
	out << quantities[coins::COIN_TYPE_DIME] << " dimes\n";
	out << quantities[coins::COIN_TYPE_QUARTER] << " quarters\n";
	out << quantities[coins::COIN_TYPE_SILVER_DOLLAR] << " silver dollars\n";
	out << "For a grand total of " << formatMoney(cents) << "!\n";
}

int main(int argc, char** argv) {
	if (argc != 2) {
		std::cout << "usage: l071 <path>\n";
		std::cout << "path: the source image for detecting coins\n";
		return -1;
	}

	std::cout << "Reading image from " << argv[1] << '\n';

	// cv: Read image from (argv[1]) in (:color) mode
	cv::Mat image = cv::imread(argv[1], cv::ImreadModes::IMREAD_COLOR);

	if (image.empty()) {
		std::cout << "No image data\n";
		return -1;
	}

	// cv: Convert the color of (image) using (:bgr2gray)
	cv::Mat gray;
	cv::cvtColor(image, gray, cv::COLOR_BGR2GRAY);

	const int GAUSSIAN_KERNEL_SIZE = 3;
	
	// cv: Blur (gray) with a kernel of size (GAUSSIAN_KERNEL_SIZE, GAUSSIAN_KERNEL_SIZE) [and std dev of (3, 3)] [using border (:default)]
	cv::medianBlur(
		gray, // input
		gray, // output
		9
	);

	const int MIN_RADIUS = 75;
	const int MAX_RADIUS = 175;
	const int MIN_DISTANCE = MIN_RADIUS;

	const int MAX_POSSIBLE_EDGE_GRADIENT = 1016;

	const int EDGE_DETECTION_THRESHOLD = MAX_POSSIBLE_EDGE_GRADIENT * 0.55 / GAUSSIAN_KERNEL_SIZE;
	const int ACCUMULATOR_THRESHOLD = 30;

	const double ACCUMULATOR_BLOCKINESS = 1;
	
	// cv: Detect edges in (gray) with thresholds (EDGE_DETECTION_THRESHOLD/2, EDGE_DETECTION_THRESHOLD)
	// Replicate the same style of edge detection used by Hough Circles
	cv::Mat edges;
	cv::Canny(gray, edges, EDGE_DETECTION_THRESHOLD / 2, EDGE_DETECTION_THRESHOLD);

	// cv::imwrite("edges.jpg", edges);

	// cv: Detect circles in [gray] at least [MIN_DISTANCE] apart with a radius from [MIN_RADIUS] to [MAX_RADIUS] and output them to [circles]
  std::vector<cv::Vec3f> circles;
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

	double estimatedPennyRadius = coins::findPennyRadius(circles, image);

	int totalValueInCents = 0;
	int coinCountsByType[6] = {0, 0, 0, 0, 0, 0};

	// for each (circle) of (circles)
	for (const auto& circle : circles) {
		// [x, y, radius] = circle
		int x = circle[0];
		int y = circle[1];
		int radius = circle[2];
		auto center = cv::Point(x, y);

		auto type = coins::estimateCoinType(radius, estimatedPennyRadius);
		// the specifications say to just use red for the circle color
		auto color = coins::COIN_COLORS[coins::COIN_TYPE_PENNY]; // coins::COIN_COLORS[type];
		// cvcolor is in BGR, coin colors are supplied in RGB
		// therefore, we make a mirror image
		auto cvcolor = cv::Scalar(color[2], color[1], color[0]);
		// cv: draw circle on (image) at (center) in color (cvcolor) [with thickness 3]
		cv::circle(image, center, radius, cvcolor, 3, cv::LINE_AA);

		totalValueInCents += coins::COIN_VALUES[type];
		coinCountsByType[type]++;
	}

	cv::imwrite("imagec.jpg", image);

	std::ofstream results("results.txt");
	writeResults(results, coinCountsByType, totalValueInCents);
	results.close();
	
	writeResults(std::cout, coinCountsByType, totalValueInCents);

	return 0;
}
