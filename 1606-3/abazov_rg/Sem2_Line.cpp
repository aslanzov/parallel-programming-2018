#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <random>

const size_t PLANE_SIZE = 10;
struct Point {
	double x;
	double y;
};

double CosVectors(double x1, double y1, double x2, double y2) {
	double dot_prod = x1 * x2 + y1 * y2;
	double x_len = sqrt(x1 * x1 + x2 * x2);
	double y_len = sqrt(y1 * y1 + y2 * y2);

	if (x_len == 0 || y_len == 0) {
		return 0;
	}
	else {
		return dot_prod / (x_len * y_len);
	}
}

unsigned JarvisStartPoint(Point* points, const size_t points_amount) {
	unsigned bottom_left = 0;

	for (size_t i = 1; i < points_amount; i++) {
		const bool below = points[i].y <= points[bottom_left].y;
		if (below || (below && points[i].x <= points[bottom_left].x)) {
			bottom_left = i;
		}
	}

	return bottom_left;
}

unsigned JarvisSecondPoint(Point* points, const size_t points_amount, const unsigned start_point) {
	unsigned second_point = start_point;
	double max_cos = -1.0;

	for (size_t i = 0; i < points_amount; i++) {
		if (i == start_point) continue;

		double cos = CosVectors(1, 0, points[i].x - points[start_point].x, points[i].y - points[start_point].y);

		if (cos >= max_cos) {
			max_cos = cos;
			second_point = i;
		}
	}

	return second_point;
}

void InitializeConvexHull(bool*& hull, const size_t size) {
	hull = new bool[size];

	for (size_t i = 0; i < size; i++) {
		hull[i] = false;
	}
}

bool* JarvisSequental( Point* points, const size_t points_amount) {
	unsigned start_point, prev_point, curr_point, next_point;
	bool *hull = nullptr;
	InitializeConvexHull(hull, points_amount);
	start_point = JarvisStartPoint(points, points_amount);
	hull[start_point] = true;

	prev_point = start_point;

	// поиск точки, имеющей наиб. полож. пол€рный угол отн. 1-й точки
	curr_point = JarvisSecondPoint(points, points_amount, start_point);
	hull[curr_point] = true;


	next_point = start_point;
	do {
		double max_cos = -1.0;

		for (size_t i = 0; i < points_amount; i++) {
			if (i != start_point && hull[i]) continue;

			double cos = CosVectors(
				points[curr_point].x - points[prev_point].x,
				points[curr_point].y - points[prev_point].y,
				points[i].x - points[curr_point].x,
				points[i].y - points[curr_point].y
			);

			if (cos >= max_cos) {
				max_cos = cos;
				next_point = i;
			}
		}

		hull[next_point] = true;
		prev_point = curr_point;
		curr_point = next_point;
	} while (next_point != start_point);
	return hull;
}

void GeneratePoints(Point*& points, const size_t points_amount) {
	points = new Point[points_amount];

	std::default_random_engine generator((unsigned)777);
	std::uniform_real_distribution<double> distribution(0, PLANE_SIZE);

	for (size_t i = 0; i < points_amount; i++) {
		points[i].x = distribution(generator);
		points[i].y = distribution(generator);
	}
}



int main(int argc, char** argv) {
	
	bool *convex_hull = nullptr;
	Point *points = nullptr;
	size_t points_amount = 0;


	struct {
		double cos;
		int rank;
	};

	std::cout << " Points number: ";
	std::cin >> points_amount;

	GeneratePoints(points, points_amount);
	convex_hull = JarvisSequental(points, points_amount);
	
	for (int i = 0; i < points_amount; i++)
	{
		if (convex_hull[i])
			std::cout << points[i].x << " " << points[i].y << std::endl;
	}


	system("pause");
	return 0;

}