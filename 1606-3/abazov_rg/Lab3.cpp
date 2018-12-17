#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <random>
#include "mpi.h"

const size_t PLANE_SIZE = 10;
const size_t POINTS_AMOUNT = 5000;


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

		double cos = CosVectors(1, 0,	points[i].x - points[start_point].x, points[i].y - points[start_point].y );

		if (cos >= max_cos) {
			max_cos = cos;
			second_point = i;
		}
	}

	return second_point;
}

void JarvisSequental(bool* hull, Point* points, const size_t points_amount) {
	unsigned start_point, prev_point, curr_point, next_point;

	//нахождение левой нижней точки
	start_point = JarvisStartPoint(points, points_amount);
	hull[start_point] = true;

	prev_point = start_point;

	// поиск точки, имеющей наим. полож. пол€рный угол отн. 1-й точки
	curr_point = JarvisSecondPoint(points, points_amount, start_point);
	hull[curr_point] = true;

	// основа
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
}

void GeneratePoints(Point*& points, const size_t points_amount) {
	points = new Point[points_amount];

	std::default_random_engine generator((unsigned)time(0));
	std::uniform_real_distribution<double> distribution(0, PLANE_SIZE);

	for (size_t i = 0; i < points_amount; i++) {
		points[i].x = distribution(generator);
		points[i].y = distribution(generator);
	}
}

void InitializeConvexHull(bool*& hull, const size_t size) {
	hull = new bool[size];

	for (size_t i = 0; i < size; i++) {
		hull[i] = false;
	}
}

bool CompareConvexHulls(bool* hull1, bool* hull2, const size_t size) {
	for (size_t i = 0; i < size; i++) {
		if (hull1[i] != hull2[i]) return false;
	}

	return true;
}

int main(int argc, char** argv) {
	int rank, proc_amount;

	bool *convex_hull_sequential = nullptr,
		*convex_hull_parallel = nullptr,
		*convex_hull_parallel_local = nullptr;
	Point *points = nullptr, *points_local = nullptr;

	size_t points_amount = 0;
	size_t points_amount_local, points_amount_local_real, points_amount_extra;

	Point start_point, prev_point, curr_point, next_point;
	int start_point_id, next_point_id, next_point_local_id;
	struct {
		double cos;
		int rank;
	} max_cos, max_cos_local;

	double time_start, time_end, time_sequential, time_parallel;

	if (argc >= 2) points_amount = atoi(argv[1]);
	if (!points_amount) points_amount = POINTS_AMOUNT;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_amount);

	if (rank == 0) {
		if (proc_amount > points_amount) {
			std::cout << "Error";
			return 1;
		}

		points_amount_extra = points_amount % proc_amount;
		if (points_amount_extra) points_amount_extra = proc_amount - points_amount_extra;
		points_amount_local = (points_amount + points_amount_extra) / proc_amount;
		points_amount_local_real = points_amount_local;

		GeneratePoints(points, points_amount + points_amount_extra);
		InitializeConvexHull(convex_hull_sequential, points_amount + points_amount_extra);
		InitializeConvexHull(convex_hull_parallel, points_amount + points_amount_extra);

		std::cout << "Amount of points: " << points_amount << std::endl;

	
		// Ћинейный
	

		time_start = MPI_Wtime();
		JarvisSequental(convex_hull_sequential, points, points_amount);
		time_end = MPI_Wtime();
		//time_sequential = time_end - time_start;

		std::cout << "Single time: " << (time_end - time_start) << std::endl;

		
		// ѕараллельный
	

		time_start = MPI_Wtime();

		start_point_id = JarvisStartPoint(points, points_amount);
		start_point = points[start_point_id];

		prev_point = start_point;
		const unsigned curr_point_id = JarvisSecondPoint(points, points_amount, start_point_id);
		curr_point = points[curr_point_id];

		convex_hull_parallel[start_point_id] = true;
		convex_hull_parallel[curr_point_id] = true;
	}

	MPI_Bcast(&points_amount_local, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&points_amount_local_real, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&points_amount_extra, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&start_point, 1, MPI_2DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
	MPI_Bcast(&start_point_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&prev_point, 1, MPI_2DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
	MPI_Bcast(&curr_point, 1, MPI_2DOUBLE_PRECISION, 0, MPI_COMM_WORLD);

	points_local = new Point[points_amount_local];
	convex_hull_parallel_local = new bool[points_amount_local];

	MPI_Scatter(points, points_amount_local, MPI_2DOUBLE_PRECISION,
		points_local, points_amount_local, MPI_2DOUBLE_PRECISION,
		0, MPI_COMM_WORLD);

	MPI_Scatter(convex_hull_parallel, points_amount_local, MPI_C_BOOL,
		convex_hull_parallel_local, points_amount_local, MPI_C_BOOL,
		0, MPI_COMM_WORLD);

	if (rank == proc_amount - 1) {
		points_amount_local_real -= points_amount_extra;
	}

	next_point_id = start_point_id;
	next_point_local_id = 0;
	do {
		max_cos_local.cos = -1;
		max_cos_local.rank = rank;

		for (size_t i = 0; i < points_amount_local_real; i++) {
			if ((i + rank * points_amount_local) != start_point_id
				&& convex_hull_parallel_local[i]) continue;

			const double cos = CosVectors(
				curr_point.x - prev_point.x,
				curr_point.y - prev_point.y,
				points_local[i].x - curr_point.x,
				points_local[i].y - curr_point.y
			);

			if (cos >= max_cos_local.cos) {
				max_cos_local.cos = cos;
				next_point_local_id = i;
			}
		}

		MPI_Allreduce(&max_cos_local, &max_cos, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
		if (rank == max_cos.rank) {
			convex_hull_parallel_local[next_point_local_id] = true;
			next_point_id = next_point_local_id + rank * points_amount_local;
			next_point = points_local[next_point_local_id];
		}
		MPI_Bcast(&next_point_id, 1, MPI_INT, max_cos.rank, MPI_COMM_WORLD);
		MPI_Bcast(&next_point, 1, MPI_2DOUBLE_PRECISION, max_cos.rank, MPI_COMM_WORLD);
		if (rank == 0) {
			convex_hull_parallel[next_point_id] = true;
		}

		prev_point = curr_point;
		curr_point = next_point;
	} while (next_point_id != start_point_id);

	if (rank == 0) {
		time_end = MPI_Wtime();
	
		std::cout << "Parall time: " << time_end - time_start << std::endl;

		const bool is_results_same = CompareConvexHulls(
			convex_hull_sequential,
			convex_hull_parallel,
			points_amount
		);

		std::cout << "Correct";
	}

	MPI_Finalize();

	return 0;
}