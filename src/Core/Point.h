// [[Rcpp::depends(BH)]]

// generic description of a point in 2D space, used for spatial hashing

#ifndef CIS_POINT_H
#define CIS_POINT_H

#include <cmath>
#include <boost/functional/hash.hpp>

template <typename T>
struct Point
{
    // coordinates of point
    T x, y;

    // less than == closer to origin
    bool operator<(const Point& other) const
    {
        return pow(x, 2) + pow(y, 2) < pow(other.x, 2) + pow(other.y, 2);
    }

    // close enough points within some tolerance
    bool operator==(const Point& other) const
    {
        return x == other.x && y == other.y;
    }

    // default definition of !=
    bool operator!=(const Point& other) const
    {
        return !(other == *this);
    }

    // default constructor
	Point() {x = 0.0; y = 0.0;}

    // constructor given specific coordinates
	Point(T inX, T inY) {x = inX; y = inY;}

    // euclidean distance between two points
	double distance(const Point& other) const
    {
		return sqrt(pow(x - other.x, 2) + pow(y - other.y, 2));
	}
};

// neccesary for boost::unordered_map
struct iequal_to : std::binary_function<Point<int>, Point<int>, bool>
{
    bool operator() (const Point<int>& p1, const Point<int>& p2) const
    {
        return p1 == p2;
    }
};

// hash function for boost::unordered_map
struct ihash : std::unary_function<Point<int>, std::size_t>
{
    std::size_t operator() (const Point<int>& p) const
    {
        // use boost to hash point coordinates
        boost::hash<int> int_hash;
        return (51 + int_hash(p.x)) * 51 + int_hash(p.y);
    }
};

#endif


