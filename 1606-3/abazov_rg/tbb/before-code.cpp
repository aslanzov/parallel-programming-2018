#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cstdlib>
#include "point.h"
#include <omp.h>
#include <iostream>

#include <random>
#include <chrono>
#include <ctime>
#include <tbb\tick_count.h>
#include <tbb\task_scheduler_init.h>

using namespace tbb;

void JarvisAlgorithm(Point* data, Point* result, int in_size, int &out_size);

int main(int argc, char* argv[]) {
	if (argc != 3) {
		std::cout << "JARVIS ALGORITHM PROGRAM\n" << "To use this program, please stick to the following pattern:\n" <<
			"solver [input] [output]" << std::endl;
		return 1;
	}
	
	int size, out_size;
	Point *in, *out;

	freopen(argv[1], "rb", stdin);

	fread(&size, sizeof(size), 1, stdin);

	in = new Point[size];
	out = new Point[size];

	fread(in, sizeof(*in), size, stdin);

	tick_count start = tick_count::now();
	tbb::task_scheduler_init init(4);
	JarvisAlgorithm(in, out, size, out_size);
	tick_count finish = tick_count::now();
	double time = (finish - start).seconds();

	freopen(argv[2], "wb", stdout);
	fwrite(&time, sizeof(time), 1, stdout);
	fwrite(&out_size, sizeof(out_size), 1, stdout);

	if (out_size > 0) {
		fwrite(out, sizeof(*out), out_size, stdout);
	}

	return 0;
}
