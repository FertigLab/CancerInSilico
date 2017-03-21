// [[Rcpp::depends(BH)]]

#ifndef CIS_LATTICE_H
#define CIS_LATTICE_H

#include <Rcpp.h>
#include <boost/unordered_map.hpp>
#include <cmath>
#include <vector>
#include <utility>

#include "Point.h"

typedef Point GridPoint;

// This data structure handles values placed in the 2D plane. It allows
// for fast iterator access to values located in a neighborhood around a point.
template <class T>
class Lattice
{
private:

    // list of values in the map (locations needed for quick deletion)
    std::vector< std::pair<GridPoint, T> > mValues;

    // map of points to index in mValueList
    boost::unordered_map<GridPoint, unsigned, ihash, iequal_to> mGrid;

    // hash a point in 2D space to the grid
    virtual GridPoint hash(const Point&) const = 0;

    // add & remove keys to the grid
    void addKey(GridPoint&, unsigned);
    void removeKey(GridPoint&);

public:

    // iterates through values within a search radius
    class local_iterator
    {    
    private:
        
        Lattice<T>* mGrid;
        GridPoint mCurrent;

        Point mCenter;
        double mRadius;

//        virtual void advance() = 0;

    public:

        local_iterator(Lattice<T>*, Point&, double, bool = false);

        local_iterator operator++(int);
        local_iterator operator++();
        bool operator!=(const local_iterator& other) const;
        T& operator*();

//       virtual void gotoEnd() = 0;
    };

    typedef typename std::vector< std::pair<GridPoint, T> >::iterator Iter;

    // full iterator
    class iterator
    {
    private:
        
        Iter it;

    public:
    
        iterator(Iter internal) : it (internal) {}

        iterator operator++() {++it; return *this;}
        iterator operator++(int) {iterator old = *this; ++it; return old;}
        bool operator!=(const iterator& other) const {return other.it != it;}
        T& operator*() {return (*it).second;}
    };

    // default constructor
    Lattice() {}

    // insert & delete values from lattice
    virtual void insert(const Point&, const T&);
    virtual void erase(const Point&);

    // update location of value
    virtual void update(const Point&, const Point&);

    // get random value from lattice
    virtual T& randomValue();

    // number of objects in lattice
    int size() const;

    // iterators
    local_iterator begin(const Point&, double) const;
    local_iterator end(const Point&, double) const;
    iterator begin() {return iterator(mValues.begin());}
    iterator end() {return iterator(mValues.end());}
};

// Adds key to the grid
// Average Case: O(1)
template <class T>
void Lattice<T>::addKey(GridPoint& key, unsigned index)
{
    if (!mGrid.insert(std::pair<GridPoint, unsigned>(key, index)).second)
    {
        throw std::invalid_argument("can't add: key already mapped\n");
    }
}

// Removes key from the grid
// Average Case: O(1)
template <class T>
void Lattice<T>::removeKey(GridPoint& key)
{
    if (mGrid.erase(key) == 0)
    {
        throw std::invalid_argument("can't remove: key is not mapped\n");
    }
}

// Inserts value in lattice
// Average Case: O(1)
template <class T>
void Lattice<T>::insert(const Point& pt, const T& val)
{
    GridPoint hashedPt = hash(pt);
    addKey(hashedPt, size());
    mValues.push_back(std::pair<GridPoint, const T>(hashedPt, val));
}

// Delete value from lattice
// Average Case: O(1)
template <class T>
void Lattice<T>::erase(const Point& pt)
{
    // delete from index grid and value list
    GridPoint hashedPt = hash(pt);
    unsigned index = mGrid.at(hashedPt);
    mGrid.erase(hashedPt);
    mValues[index] = mValues.back();
    mValues.pop_back();

    // update key of moved object in previous step
    mGrid.erase(mValues[index].first);
    mGrid.insert(std::pair<GridPoint, unsigned>(mValues[index].first, index));
}

// Update a value after it moves
// Average Case: O(1)
template <class T>
void Lattice<T>::update(const Point& oldPt, const Point& newPt)
{
    T val = mValues[mGrid.at(hash(oldPt))].second;
    
    erase(oldPt);
    insert(newPt, val);
}

// Return a random value in the lattice
// O(1)
template <class T>
T& Lattice<T>::randomValue()
{
    return mValues[floor(R::runif(0, size()))].second;
}

// return size of lattice
template <class T>
int Lattice<T>::size() const
{
    return mValues.size();
}

#endif





