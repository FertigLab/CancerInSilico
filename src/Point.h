// [[Rcpp::depends(BH)]]

#ifndef POINT_HPP
#define POINT_HPP

#include <cmath>
#include <boost/functional/hash.hpp>

typedef struct point {

    double x;
    double y;

    bool operator<(const struct point& other) const {
        return pow(x, 2) + pow(y, 2) < pow(other.x, 2) + pow(other.y, 2);
    }

    bool operator==(const struct point& other) const {
        return (x < other.x + 0.0001 && x > other.x - 0.0001)
				&& (y < other.y + 0.0001 && y > other.y - 0.0001);
    }

    bool operator!=(const struct point& other) const {
        return !(other == *this);
    }

	point() { x = 0.0; y = 0.0;}

	point(double in_x, double in_y) {
		x = in_x;
		y = in_y;
	}

	double dist(const struct point& other) const {
		return pow(pow(x - other.x,2) + pow(y - other.y,2),0.5);
	}

	void operator=(const struct point& other) {
		x = other.x;
		y = other.y;
	}

} Point;

struct iequal_to
    : std::binary_function<Point, Point, bool> {

    bool operator()(const Point &p1, const Point &p2) const {
        return p1.x == p2.x && p1.y == p2.y;
    }

};

struct ihash
    : std::unary_function<Point, std::size_t> {

    std::size_t operator()(const Point &p) const {
        boost::hash<int> int_hash;
        return (51 + int_hash(p.x)) * 51 + int_hash(p.y);
    }

};

#endif
