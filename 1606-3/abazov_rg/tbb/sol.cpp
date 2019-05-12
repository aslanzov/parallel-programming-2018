#define _CRT_SECURE_NO_WARNINGS

#include <cmath>
#include "point.h"
#include <tbb\task_scheduler_init.h>
#include <tbb\tbb.h>

using namespace tbb;

double Cos(double x1, double y1, double x2, double y2);

class CosCompare {
public:
	explicit CosCompare(Point *data, int prev_p, int cur_p, bool *is_used) :
		_data(data), _prev_p(prev_p), _cur_p(cur_p), _is_used(is_used), _min_cos(1), _next_p(0) {}

	CosCompare(const CosCompare &c, split) :
		_data(c._data), _prev_p(c._prev_p), _cur_p(c._cur_p), _is_used(c._is_used), _min_cos(1), _next_p(0) {}

	void operator() (const blocked_range<int> &r) {
		for (int i = r.begin(); i != r.end(); i++) {
			if (!_is_used[i]) {
				double cur_cos = Cos(_data[_prev_p].x - _data[_cur_p].x, _data[_prev_p].y - _data[_cur_p].y,
					_data[i].x - _data[_cur_p].x, _data[i].y - _data[_cur_p].y);

				if (cur_cos < _min_cos) {
					_next_p = i;
					_min_cos = cur_cos;
				}
			}
		}
	}

	void join(const CosCompare &c) {
		if (c._min_cos < _min_cos)
			_next_p = c._next_p;
	}

	int Result() {
		return _next_p;
	}
private:
	int _prev_p, _cur_p, _next_p;
	double _min_cos;
	Point *_data;
	bool* _is_used;
};

double Cos(double x1, double y1, double x2, double y2) {
	double scalar_prod = x1 * x2 + y1 * y2;
	double a1_len = sqrt(x1 * x1 + y1 * y1);
	double a2_len = sqrt(x2 * x2 + y2 * y2);

	double mult_of_len = a1_len * a2_len;

	if (mult_of_len == 0)
		return 0;
	else
		return scalar_prod / mult_of_len;
}

void JarvisAlgorithm(Point* data, Point* result, int in_size, int &out_size) {
	int first_p, prev_p, cur_p, next_p;
	bool* is_used;
	int *indexes;
	int k = 0;
	double max_cos;

	if (in_size < 3) {
		out_size = 0;
		return;
	}

	indexes = new int[in_size];
	is_used = new bool[in_size];

	first_p = 0;

	for (int i = 0; i < in_size; i++) {
		indexes[i] = -1;
		is_used[i] = false;
	}

	//Looking for a first point
	for (int i = 0; i < in_size; i++) {
		if ((data[i].y < data[first_p].y) ||
			(data[i].y == data[first_p].y) && (data[i].x < data[first_p].x))
			first_p = i;
	}

	indexes[k++] = first_p;

	//Looking for a second point
	max_cos = -1;

	for (int i = 0; i < in_size; i++) {
		if (!is_used[i]) {
			double cur_cos = Cos(1, 0, data[i].x - data[first_p].x,
				data[i].y - data[first_p].y);

			if (cur_cos > max_cos) {
				cur_p = i;
				max_cos = cur_cos;
			}
		}
	}

	indexes[k++] = cur_p;
	is_used[cur_p] = true;

	//The remaining steps of computing the convex hull
	prev_p = first_p;

	task_scheduler_init init(6);

	do {
		CosCompare main_compare (data, prev_p, cur_p, is_used);
		parallel_reduce(blocked_range<int>(0, in_size, 10), main_compare);
		next_p = main_compare.Result();

		if (next_p != first_p) {
			indexes[k++] = next_p;
			is_used[next_p] = true;
			prev_p = cur_p;
			cur_p = next_p;
		}
	} while (next_p != first_p);

	init.terminate();

	out_size = k;

	for (int i = 0; i < k; i++) {
		result[i] = data[indexes[i]];
	}
}

