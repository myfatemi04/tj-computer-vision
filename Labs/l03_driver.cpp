#ifndef LAB_PART
#define LAB_PART -1
#endif

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

	// array of booleans: whether or not to test out each part
	// if parts[2] is true, it will test the lab 3 part 2 code.
	// you can test out multiple parts.
	bool parts[] = {false, false, false, false};

	// how many cycles to run when timing the parts
	int cycles = 10;

	// if n > 0, generates n points and stores them in a file
	int npointsToGenerate = 0;

	// the file to store the generated points to and read from
	const char* file = "points.txt";

	// if no arguments are present, print out a help message
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

	// read options from the args
	for (int i = 1; i < argc; i++) {
		// 'g' to generate N points
		if (argv[i][0] == 'g') {
			npointsToGenerate = atoi(argv[i] + 1);

		// 'p' to specify which parts to run
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

		// 'c' to specify how many cycles to run
		} else if (argv[i][0] == 'c') {
			cycles = atoi(argv[i] + 1);

		// 'f' to specify the input file
		} else if (argv[i][0] == 'f') {
			file = argv[i] + 1;
		}
	}

	if (npointsToGenerate > 0) {
		std::vector<Point> points = generatePoints(npointsToGenerate);
		savePoints(points, file);
		std::cout << "Generated " << npointsToGenerate << " points" << std::endl;
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
