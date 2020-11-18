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
	int npointsToGenerate = 0; // if 0, don't generate any
	const char* file = "points.txt";

	if (argc == 1) {
		const char* help = 
"Usage: l03_driver gN p123 cN iFile\n"
"    gN: Generate N points\n"
"    p123: Time part 1, 2, and/or 3. Can do p1, p23, etc.\n"
"    cN: Cycle N times for the timer. Defaults to 10 cycles.\n"
"    fFile: Use the file File for generating/reading points.\n"
"           Defaults to \"points.txt\".\n";
		std::cout << help;
		return 0;
	}

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == 'g') {
			npointsToGenerate = atoi(argv[i] + 1);
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
		} else if (argv[i][0] == 'f') {
			file = argv[i] + 1;
		}
	}

	if (npointsToGenerate > 0) {
		std::vector<Point> points = generatePoints(npointsToGenerate);
		savePoints(points, file);
	}

  auto points = readPoints(file);

  std::ofstream outfile("results.txt");
	std::cout << std::setprecision(17);
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
