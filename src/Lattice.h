// [[Rcpp::depends(BH)]]

#ifndef CIS_LATTICE_H
#define CIS_LATTICE_H

#include <Rcpp.h>
#include <boost/unordered_map.hpp>
#include <cmath>
#include <vector>
#include <utility>
#include <exception>

#include "Point.h"

typedef Point<int> GridPoint;
template <class T> class BaseLocalIterator;

// This data structure handles values placed in the 2D plane. It allows for
// fast iterator access to values located in a neighborhood around a point.
template <class T>
class Lattice
{

friend class BaseLocalIterator<T>;

protected:

    // list of values in the map (locations needed for quick deletion)
    std::vector< std::pair<GridPoint, T> > mValues;

    // map of points to index in mValueList
    boost::unordered_map<GridPoint, unsigned, ihash, iequal_to> mGrid;

    // hash a point in 2D space to the grid
    virtual GridPoint hash(const Point<double>&) const = 0;

    // add & remove keys to the grid
    void addKey(GridPoint&, unsigned);
    void removeKey(GridPoint&);

public:

    /***************** standard iterator ******************/
    typedef typename std::vector< std::pair<GridPoint, T> >::iterator Iter;

    class iterator
    {
    private:
        
        Iter it;

    public:
    
        iterator(Iter internal) : it(internal) {}
        iterator(const iterator& other) {it = other.it;}
        ~iterator() {}

        iterator& operator=(const iterator& itr) {it = itr.it; return *this;}
        iterator& operator++() {++it; return *this;}
        iterator operator++(int) {iterator tmp(*this); ++it; return tmp;}
        bool operator!=(const iterator& itr) const {return itr.it != it;}
        T& operator*() {return (*it).second;}
    };
    /******************************************************/
  
    /****************** local iterator ********************/
    class local_iterator
    {
    protected:
    
        BaseLocalIterator<T>* baseIter;

    public:
        
        local_iterator(BaseLocalIterator<T>* b)
            : baseIter(b) {}
        local_iterator(const local_iterator& i)
            {baseIter = i.baseIter->newCopy();}
        ~local_iterator()
            {delete baseIter;}
        
        local_iterator& operator=(const local_iterator& i)
            {baseIter = i.baseIter; return *this;}
        local_iterator& operator++()
            {baseIter->operator++(); return *this;}
        bool operator!=(const local_iterator& i) const
            {return baseIter->operator!=(i.baseIter);}
        T& operator*()
            {return baseIter->operator*();}
    };
    /******************************************************/

    // default constructor
    Lattice() {}

    // insert & delete values from lattice
    virtual void insert(const Point<double>&, const T&);
    virtual void erase(const Point<double>&);

    // update location of value
    virtual void update(const Point<double>&, const Point<double>&);

    // get random value from lattice
    virtual T& randomValue();

    // number of objects in lattice
    int size() const;

    // iterators
    virtual local_iterator lbegin(const Point<double>&, double) = 0;
    virtual local_iterator lend(const Point<double>&, double) = 0;
    iterator begin() {return iterator(mValues.begin());}
    iterator end() {return iterator(mValues.end());}
};

template <class T>
class BaseLocalIterator
{    
protected:
    
    Lattice<T>* mLattice;
    GridPoint mCurrent;

    Point<double> mCenter;
    double mRadius;

public:

    BaseLocalIterator(Lattice<T>* l, const Point<double>& c, double r)
        : mLattice(l), mCenter(c), mRadius(r) {}

    virtual BaseLocalIterator* newCopy() const = 0;
    virtual void operator++() = 0;

    T& operator*()
        {return mLattice->mValues[mLattice->mGrid.at(mCurrent)];}

    bool operator!=(const BaseLocalIterator* i) const
    {
        if ((mCenter != i->mCenter) || (mRadius != i->mRadius))
        {
            throw std::invalid_argument("invalid iterator comparison");
        }
        else
        {
            return mCurrent != i->mCurrent;
        }
    }        
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
void Lattice<T>::insert(const Point<double>& pt, const T& val)
{
    GridPoint hashedPt = hash(pt);
    addKey(hashedPt, size());
    mValues.push_back(std::pair<GridPoint, const T>(hashedPt, val));
}

// Delete value from lattice
// Average Case: O(1)
template <class T>
void Lattice<T>::erase(const Point<double>& pt)
{
    // delete from index grid and value list
    GridPoint hashedPt = hash(pt);
    unsigned index = mGrid.at(hashedPt);
    mGrid.erase(hashedPt);

    // update key of moved object
    if (index < mValues.size() - 1)
    {
        mGrid.erase(mValues.back().first);
        mGrid.insert(std::pair<GridPoint, unsigned>(mValues.back().first,
            index));
    }
    
    // delete value
    mValues[index] = mValues.back();
    mValues.pop_back();
}

// Update a value after it moves
// Average Case: O(1)
template <class T>
void Lattice<T>::update(const Point<double>& oldPt, const Point<double>& newPt)
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





