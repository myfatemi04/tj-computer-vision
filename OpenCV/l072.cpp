#include <opencv2/core.hpp>

typedef struct {
	double x, y;
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

// Returns an array of eight vertices.
Location3 *getVerticesOfCube(const Cube& cube) {
	return nullptr;
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
	// Get the vector from the viewer to the point
	auto relativeVector = getVectorFromAToB(perspective.location, point);
	
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
	// components of the relative vector. We can ignore the Z component.
	auto projectionPlaneDistance = camera.projectionPlaneDistance;
	auto u = vectorFromPerspective.x * projectionPlaneDistance;
	auto v = vectorFromPerspective.y * projectionPlaneDistance;
	
	Location2 projectedLocation {
		u, v
	};

	return projectedLocation;
}

void renderCube(const cv::Mat& out, const Cube& cube) {
	auto vertices = getVerticesOfCube(cube);
	
}

int main() {

}