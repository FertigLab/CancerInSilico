// [[Rcpp::depends(BH)]]

#ifndef CIS_POINT_H
#define CIS_POINT_H

#include <cmath>
#include <boost/functional/hash.hpp>

struct Point
{
    /* coordinates of point */
    double x, y;

    /* "less than" means closer to origin */
    bool operator<(const Point& other) const
    {
        return pow(x, 2) + pow(y, 2) < pow(other.x, 2) + pow(other.y, 2);
    }

    /* close enough points within some tolerance */
    bool operator==(const Point& other) const
    {
        return (x < other.x + 0.0001 && x > other.x - 0.0001)
            && (y < other.y + 0.0001 && y > other.y - 0.0001);
    }

    /* default definition of != */
    bool operator!=(const Point& other) const
    {
        return !(other == *this);
    }

    /* default constructor */
	Point() { x = 0.0; y = 0.0;}

    /* constructor given specific coordinates */
	Point(double inX, double inY)
    {
		x = inX;
		y = inY;
	}

    /* euclidean distance between two points */
	double distance(const Point& other) const
    {
		return pow(pow(x - other.x,2) + pow(y - other.y,2), 0.5);
	}
};

/* neccesary for boost::unordered_map */
struct iequal_to : std::binary_function<Point, Point, bool>
{
    bool operator() (const Point& p1, const Point& p2) const
    {
        return p1 == p2;
    }
};

/* hash function - used in SpatialHash */
struct ihash : std::unary_function<Point, std::size_t>
{
    std::size_t operator() (const Point& p) const
    {
        /* use boost to hash point coordinates */
        boost::hash<int> int_hash;
        return (51 + int_hash(p.x)) * 51 + int_hash(p.y);
    }
};

#endif


