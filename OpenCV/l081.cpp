#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#define SHOW_IMAGE_DEBUG_MESSAGES false

#if SHOW_IMAGE_DEBUG_MESSAGES
#  define idbg(x) std::cout << x
#else
#  define idbg(x)
#endif

class Location2 {
	private:
		double x, y;

	public:
		Location2(double x, double y): x(x), y(y) {}

		double getX() const {
			return x;
		}
		
		double getY() const {
			return y;
		}

		Location2 scale(double amount) const {
			return Location2 { x * amount, y * amount };
		}

		std::string toPointString() const {
			std::stringstream s;
			s << '(' << x << ", " << y << ')';
			return s.str();
		}

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
};

typedef double Orientation2;

typedef struct {
	Location2 location;
	Orientation2 orientation;
} Position2;

class Location3 {
	private:
		double x, y, z;

	public:
		Location3(): x(0), y(0), z(0) {}

		Location3(double x, double y, double z): x(x), y(y), z(z) {}

		double getX() const {
			return x;
		}

		double getY() const {
			return y;
		}

		double getZ() const {
			return z;
		}

		void rotateAroundZ(double radians) {
			double _x = x;
			double _y = y;
			x = cos(radians) * _x - sin(radians) * _y;
			y = sin(radians) * _x + cos(radians) * _y;
		}

		void translate(const Location3& amount) {
			x += amount.x;
			y += amount.y;
			z += amount.z;
		}

		std::string toVectorString() const {
			std::stringstream s;
			s << '<' << x << ", " << y << ", " << z << '>';
			return s.str();
		}
		
		std::string toPointString() const {
			std::stringstream s;
			s << '(' << x << ", " << y << ", " << z << ')';
			return s.str();
		}

		Location3 vectorTo(const Location3& other) const {
			return Location3 {
				other.x - x,
				other.y - y,
				other.z - z
			};
		}

		Location3 rotatedAroundZ(double rotation) const {
			double cosTheta = cos(rotation);
			double sinTheta = sin(rotation);

			Location3 rotated {
				cosTheta * x - sinTheta * y,
				sinTheta * x + cosTheta * y,
				z
			};

			return rotated;
		}
};

typedef double Orientation3;

typedef struct {
	Location3 location;
	Orientation3 orientation;
} Position3;

class Cube {
	private:
		Position3 position;
		double sideLength;

	public:
		Cube(Position3 position, double sideLength): position(position), sideLength(sideLength) {}
		
		Position3 &getPosition() {
			return position;
		}

		Location3 *getVertices() const {
			double scaleFactor = sideLength / 2;
			double rotation = position.orientation;

			auto locations = new Location3[8];
			for (int i = 0; i < 8; i++) {
				// Points are ordered by YXZ in bits
				int y = (i >> 2) & 1;
				int x = (i >> 1) & 1;
				int z = (i >> 0) & 1;

				auto vertex = Location3 {
					(x * sideLength) - scaleFactor,
					(y * sideLength) - scaleFactor,
					(z * sideLength) - scaleFactor
				};

				vertex.rotateAroundZ(rotation);

				vertex.translate(position.location);

				locations[i] = vertex;
			}

			return locations;
		}
};

typedef struct {
	Position3 position;
	double projectionPlaneDistance;
} Camera3;

typedef int VertexID;

typedef struct {
	VertexID a;
	VertexID b;
} Edge;

Location3 getVectorFromPerspective(const Location3& point, const Position3& perspective) {
	// Get the vector from the viewer to the point
	auto relativeVector = perspective.location.vectorTo(point);
	// Then, rotate it by the viewer's rotation
	auto zRotation = perspective.orientation;
	auto rotatedRelativeVector = relativeVector.rotatedAroundZ(-zRotation);

	return rotatedRelativeVector;
}

Location2 getProjectedLocation(const Location3& point, const Camera3& camera) {
	// Get the vector from the camera to the point we want to project
	auto vectorFromPerspective = getVectorFromPerspective(point, camera.position);
	// Take the Y and Z values as U and V.
	// We simply multiply the distance to the projection plane by the Y and Z
	// components of the relative vector, scaled so that X = 1.
	auto vectorFromPerspectiveScaled = Location2 {
			vectorFromPerspective.getY(), // y is side to side
			vectorFromPerspective.getZ()  // z is up and down
	}.scale(1 / vectorFromPerspective.getX()); // x is out from the screen

	auto projectionPlaneDistance = camera.projectionPlaneDistance;

	auto projectedLocation = vectorFromPerspectiveScaled.scale(projectionPlaneDistance);

	idbg("Projected location for " << point.toPointString() << " (perspective " << vectorFromPerspective.toVectorString() << ')');
	idbg(" is at " << projectedLocation.toPointString() << '\n');

	return projectedLocation;
}

void renderCube(cv::Mat& out, const Camera3& camera, const Cube& cube) {
	auto vertices = cube.getVertices();
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
	auto red = cv::Scalar(255, 0, 0);

	for (int i = 0; i < 12; i++) {
		auto edge = edges[i];
		auto firstVertexID = edge.a;
		auto secondVertexID = edge.b;

		auto firstVertexProjectedLocation = projectedVertices[firstVertexID];
		auto secondVertexProjectedLocation = projectedVertices[secondVertexID];

		auto firstVertexCentered = firstVertexProjectedLocation.toCenteredPoint(out.size);
		auto secondVertexCentered = secondVertexProjectedLocation.toCenteredPoint(out.size);

		// Draw line from first vertex to second vertex
		cv::line(out, firstVertexCentered, secondVertexCentered, i == 0 ? red : white);
	}
}

int main() {
	const int IMAGE_WIDTH = 1200;
	const int IMAGE_HEIGHT = 800;

	auto imageSize = cv::Size2i(IMAGE_WIDTH, IMAGE_HEIGHT);
	auto cube = Cube(Position3 { Location3(0, 0, 0), 1 }, 1);
	auto camera = Camera3 { Position3 { Location3 { -5, 0, 0 }, 0 }, 3 };
	auto image = cv::Mat(imageSize, CV_8UC3);

	cv::VideoWriter writer;

	writer.open("rotation.avi", cv::VideoWriter::fourcc('M', 'J', 'P', 'G'), 30, imageSize, true);

	for (int i = 0; i < 30 * 3.14159265 * 8; i++) {
		image = cv::Scalar(0, 0, 0);

		renderCube(image, camera, cube);

		writer.write(image);

		cube.getPosition().orientation = (i / (30 * 3.14159265));
	}

	writer.release();
}