#define LAB_03_MAIN

#include "l03core.cpp"
#include "l031.cpp"
#include "l032.cpp"
#include "l033.cpp"

	/**
 * Main method. Sets the random seed to the current time.
 * Also takes in an argument for the number of points to use.
 * If not provided, the default number of points is 10000.
 * Reads the list of points from the default file.
 */
int main(int argc, const char* argv[]) {
  std::srand(time(NULL));

	bool parts[] = {false, false, false, false};
	int cycles = 10;
	const char* file = "points.txt";

	if (argc == 1) {
		std::cout << "Usage: l03_driver gN p123 cN iFile\n";
		std::cout << "    gN: Generate N points\n";
		std::cout << "    p123: Time part 1, 2, and/or 3. Can do p1, p23, etc.\n";
		std::cout << "    cN: Cycle N times for the timer.\n";
		std::cout << "    iFile: Use the file File for generating/reading points.\n";
		std::cout << "           Defaults to \"points.txt\".\n";
		return 0;
	}

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == 'g') {
			int npoints = atoi(argv[i] + 1);
			std::vector<Point> points = generatePoints(npoints);
			savePoints(points);
		} else if (argv[i][0] == 'p') {
			const char* p = argv[i];
			while (*p != '\0') {
				if (*p == '1') {
					parts[1] = true;
				} else if (*p == '2') {
					parts[2] = true;
				} else if (*p == '3') {
					parts[3] = true;
				}
				p++;
			}
		} else if (argv[i][0] == 'c') {
			cycles = atoi(argv[i] + 1);
		} else if (argv[i][0] == 'i') {
			file = argv[i] + 1;
		}
	}

  auto points = readPoints();

  std::ofstream outfile("results.txt");
	std::cout << "Found " << points.size() << " points\n";
	outfile << "Found " << points.size() << " points\n";
	if (parts[1]) {
		timer(points, outfile, part1, "Brute Force", cycles);
	}
	if (parts[2]) {
		timer(points, outfile, part2, "Recursive", cycles);
	}
	if (parts[3]) {
		timer(points, outfile, part3, "Recursive Optimized", cycles);
	}
  outfile.close();
}
