// [[Rcpp::depends(BH)]]

#ifndef POINT_HPP
#define POINT_HPP

#include <cmath>
#include <boost/functional/hash.hpp>

struct Point
{
    /* coordinates of point */
    double x, y;

    /* "less than" means closer to origin */
    bool operator<(const struct point& other) const
    {
        return pow(x, 2) + pow(y, 2) < pow(other.x, 2) + pow(other.y, 2);
    }

    /* close enough points within some tolerance */
    bool operator==(const struct point& other) const
    {
        return (x < other.x + 0.0001 && x > other.x - 0.0001)
            && (y < other.y + 0.0001 && y > other.y - 0.0001);
    }

    /* default definition of != */
    bool operator!=(const struct point& other) const
    {
        return !(other == *this);
    }

    /* default constructor */
	Point() { x = 0.0; y = 0.0;}

    /* constructor given specific coordinates */
	point(double in_x, double in_y)
    {
		x = in_x;
		y = in_y;
	}

    /* euclidean distance between two points */
	double distance(const struct point& other) const
    {
		return pow(pow(x - other.x,2) + pow(y - other.y,2), 0.5);
	}
};

/*struct iequal_to
    : std::binary_function<Point, Point, bool> {

    bool operator()(const Point &p1, const Point &p2) const {
        return p1.x == p2.x && p1.y == p2.y;
    }

};*/

/* hash function - used in SpatialHash */
struct ihash : std::unary_function<Point, std::size_t>
{
    std::size_t operator()(const Point &p) const
    {
        /* use boost to hash point coordinates */
        boost::hash<int> int_hash;
        return (51 + int_hash(p.x)) * 51 + int_hash(p.y);
    }
};

#endif


