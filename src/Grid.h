// [[Rcpp::depends(BH)]]

#ifndef CIS_SPATIAL_HASH_HPP
#define CIS_SPATIAL_HASH_HPP

#include <R.h>
#include <Rcpp.h>
#include <boost/unordered_map.hpp>
#include <cmath>
#include <vector>
#include <utility>

#include "Point.h"

typedef Point GridPoint;

/*
    This class handles objects placed in the 2-D plane. It allows
    for fast iterator access to objects located in a neighborhood 
    around a point. 
*/
template <class T>
class Grid
{

private:

    /* map of grid point - value (note ihash is the internal hash) */
    boost::unordered_map<GridPoint, T*, ihash, iequal_to> mValueGrid;

    /* list of values in the map - used to get random value in O(1) */
    std::vector<T*> mValueList;

    /* map of points to index in mValueList - makes deletion O(1) */
    boost::unordered_map<GridPoint, unsigned, ihash, iequal_to> mIndexGrid;

public:

    class circular_iterator;
    typedef typename std::vector<T*>::iterator iterator.

    /* construct with fixed width between grid lines */
    Grid(double);

    /* insert & delete objects from grid */
    void insert(const Point&, const T*);
    void erase(const Point&, const T*);

    /* update location of object */
    void update(const Point&, const Point&);

    /* get random object from map */
    T* randomValue() const;

    /* number of objects in map */
    int size() const;

    /* iterators */
    circular_iterator begin(const Point&, double) const;
    circular_iterator end(const Point&, double) const;
    iterator begin() const;
    iterator end() const;
};

template <class T>
class Grid::circular_iterator
{
private:
    
    struct Box {double left, right, top, bottom;} mSearchRegion;

    Grid<T>& mGrid;
    GridPoint mCurrent;

    Point mCenter;
    double mRadius;

    void constructRegion(double radius);
    void advance();

public:

    circular_iterator(Grid<T>&, Point&, double, bool = false);

    circular_iterator operator++(int);
    circular_iterator operator++();
    T& operator*();
    T* operator&();
    GridPoint location();

    bool operator!=(const circular_iterator& other) const;
    void gotoEnd();
};

#include "Grid.cpp"

#endif
