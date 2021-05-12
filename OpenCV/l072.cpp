#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

typedef struct {
	double x, y;

	cv::Point toPoint() const {
		return cv::Point(x, y);
	}

	cv::Point toCenteredPoint(const cv::MatSize& size) const {
		int width = size[1];
		int height = size[0];
		int centerX = width / 2;
		int centerY = height / 2;
		int projectionToImageScalar = (height / 4);
		return cv::Point((x * projectionToImageScalar) + centerX, (y * projectionToImageScalar) + centerY);
	}
} Location2;

typedef double Orientation2;

typedef struct {
	Location2 location;
	Orientation2 orientation;
} Position2;

typedef struct {
	double x, y, z;
} Location3;

typedef double Orientation3;

typedef struct {
	Location3 location;
	Orientation3 orientation;
} Position3;

typedef struct {
	Position3 position;
	double sideLength;
} Cube;

typedef struct {
	Position3 position;
	double projectionPlaneDistance;
} Camera3;

typedef int VertexID;

typedef struct {
	VertexID a;
	VertexID b;
} Edge;

// Returns an array of eight vertices.
Location3 *getVerticesOfCube(const Cube& cube) {
	double scaleFactor = cube.sideLength / 2;
	double rotation = cube.position.orientation;

	auto locations = new Location3[8];
	for (int i = 0; i < 8; i++) {
		// Points are ordered by YZX in bits
		int x = (i >> 0) & 1;
		int z = (i >> 1) & 1;
		int y = (i >> 2) & 1;

		locations[i] = Location3 {
			(x * cube.sideLength) - scaleFactor,
			(y * cube.sideLength) - scaleFactor,
			(z * cube.sideLength) - scaleFactor
		};
	}

	return locations;
}

Location3 getVectorFromAToB(const Location3& b, const Location3& a) {
	return {
		b.x - a.x,
		b.y - a.y,
		b.z - a.z
	};
}

// Rotates the vector around the Z axis.
Location3 getVectorRotatedAroundZ(const Location3& a, double rotation) {
	double cosTheta = cos(rotation);
	double sinTheta = sin(rotation);

	Location3 rotated {
		cosTheta * a.x - sinTheta * a.y,
		sinTheta * a.x + cosTheta * a.y,
		a.z
	};

	return rotated;
}

Location3 getVectorFromPerspective(const Location3& point, const Position3& perspective) {
	std::cout << "Getting vector from location to point\n";
	// Get the vector from the viewer to the point
	auto relativeVector = getVectorFromAToB(perspective.location, point);

	std::cout << "Relative vector: ";
	std::cout << "<" << relativeVector.x << ", " << relativeVector.y << ", " << relativeVector.z << ">";
	std::cout << '\n';
	
	std::cout << "Rotating vector around Z\n";
	// Then, rotate it by the viewer's rotation
	auto zRotation = perspective.orientation;
	auto rotatedRelativeVector = getVectorRotatedAroundZ(relativeVector, -zRotation);

	return rotatedRelativeVector;
}

Location2 getProjectedLocation(const Location3& point, const Camera3& camera) {
	// Get the vector from the camera to the point we want to project
	auto vectorFromPerspective = getVectorFromPerspective(point, camera.position);
	// Take the X and Z values as U and V.
	// We simply multiply the distance to the projection plane by the X and Y
	// components of the relative vector, scaled so that Z = -1.
	auto vectorFromPerspectiveScaled = Location2 {
			vectorFromPerspective.x / vectorFromPerspective.z,
			vectorFromPerspective.y / vectorFromPerspective.z
	};
	auto projectionPlaneDistance = camera.projectionPlaneDistance;
	auto u = vectorFromPerspectiveScaled.x * projectionPlaneDistance;
	auto v = vectorFromPerspectiveScaled.y * projectionPlaneDistance;
	
	Location2 projectedLocation {
		u, v
	};

	std::cout << "Projected location for ";
	std::cout << "(" << point.x << ", " << point.y << ", " << point.z << ")";
	std::cout << " is at (" << u << ", " << v << ")\n";

	return projectedLocation;
}

void renderCube(cv::Mat& out, const Camera3& camera, const Cube& cube) {
	std::cout << "Getting vertices of cube\n";
	auto vertices = getVerticesOfCube(cube);
	auto edges = new Edge[12] {
		{0, 1},
		{1, 3},
		{3, 2},
		{2, 0},

		{4, 5},
		{5, 7},
		{7, 6},
		{6, 4},

		{0, 4},
		{1, 5},
		{2, 6},
		{3, 7}
	};

	auto projectedVertices = new Location2[8] {
		getProjectedLocation(vertices[0], camera),
		getProjectedLocation(vertices[1], camera),
		getProjectedLocation(vertices[2], camera),
		getProjectedLocation(vertices[3], camera),
		getProjectedLocation(vertices[4], camera),
		getProjectedLocation(vertices[5], camera),
		getProjectedLocation(vertices[6], camera),
		getProjectedLocation(vertices[7], camera),
	};

	auto white = cv::Scalar(255, 255, 255);

	for (int i = 0; i < 12; i++) {
		auto edge = edges[i];
		auto firstVertexID = edge.a;
		auto secondVertexID = edge.b;

		auto firstVertexProjectedLocation = projectedVertices[firstVertexID];
		auto secondVertexProjectedLocation = projectedVertices[secondVertexID];

		auto firstVertexCentered = firstVertexProjectedLocation.toCenteredPoint(out.size);
		auto secondVertexCentered = secondVertexProjectedLocation.toCenteredPoint(out.size);

		// Draw line from first vertex to second vertex
		cv::line(out, firstVertexCentered, secondVertexCentered, white);
	}
}

int main() {
	const int IMAGE_WIDTH = 1200;
	const int IMAGE_HEIGHT = 800;
	auto image = cv::Mat(cv::Size2i(IMAGE_WIDTH, IMAGE_HEIGHT), CV_8UC1);
	auto cube = Cube { Position3 { Location3 { 0, 0, 0 }, 0 }, 1 };
	auto camera = Camera3 { Position3 { Location3 { 0, 0, -1 }, 0 }, 1 };

	std::cout << "Rendering cube\n";

	renderCube(image, camera, cube);

	cv::imwrite("cube.png", image);
}