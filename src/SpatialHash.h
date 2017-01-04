// [[Rcpp::depends(BH)]]

#ifndef SPATIAL_HASH_HPP
#define SPATIAL_HASH_HPP

#include <R.h>
#include <Rcpp.h>
#include <boost/unordered_map.hpp>
#include <cmath>
#include <vector>

#include "Point.h"

/*
    This class handles objects placed in the 2-D plane. It allows
    for fast iterator access to objects located in a neighborhood 
    around a point. First a grid is established based on the minimum
    radius of an object of type T. The grid ensures each square can 
    have at most one occupant. This allows for quick hashing between
    objects and grid points. The grid points are labeled by indices, 
    not actual coordinates, i.e. (0,0) - (1,0) - (2,0) - etc..
*/
template <class T>
class SpatialHash
{

friend class TestSpatialHash;
friend class TestCellPopulation; //bad

private:

    /* map of grid point - value (note ihash is the internal hash) */
    boost::unordered_map<Point, T*, ihash, iequal_to> mHashMap;

    /* list of values in the map - used to get random value in O(1) */
    std::vector<T*> mValueList;

    /* map of points to index in mValueList - makes deletion O(1) */
    boost::unordered_map<Point, unsigned, ihash, iequal_to> mIndexMap;

    /* space between grid lines */
    double mBucketSize;

    /* get nearest coordinate on grid */
    Point hash(Point);

    /* add & remove keys from the hash map */
    void removeKey(Point);
    void addKey(Point, T*);

public:

    /* iterator around a point */
    class circular_iterator
    {

    private:
    
        /* square boundaries of search region */
        struct box
        {
            double left, right, top, bottom;
        } mSearchRegion;

        /* reference to hash map */
        SpatialHash<T>& mHash;

        /* current location of iterator */
        Point mCurrent;

        /* needed for establishing equality */
        Point mCenter;
        double mRadius;

        /* create search region */
        void constructRegion(double radius)
        {
            /* construct region */
            mSearchRegion.left = mHash.hash(mCenter.x - mRadius);
            mSearchRegion.right = mHash.hash(mCenter.x + mRadius);
            mSearchRegion.bottom = mHash.hash(mCenter.y - mRadius);
            mSearchRegion.top = mHash.hash(mCenter.y + mRadius);
        }

        /* advance to next location */
        void advance()
        {
            do
            {
                /* if outside region, stop iterator */
                if (mCurrent.x > mSearchRegion.right)
                {
                    break;
                }

                /* move down one */
                mCurrent.y -= 1.0;
                if (mCurrent.y < mSearchRegion.bottom)
                {
                    /* if outside region, move to top, over one */
                    mCurrent.y = mSearchRegion.top;
                    mCurrent.x += 1.0;
                }
            /* iterate through region until point (not center) is found */
            } while (mHash->mHashMap.count(mCurrent) == 0 
                || mCurrent == mCenter);
        }

    public:

        /* constructor - hash, point, and radius - bool for end iterator */
        circular_iterator(SpatialHash<T>& hash, Point& center,
            double radius, bool end = false)
        {
            /* store variables */
            mHash = hash;   
            mCenter = center;
            mRadius = radius;

            /* calculate search region */
            constructRegion();
              
            /* get starting position of iterator */
            mCurrent.x = mSearchRegion.left - 1.0;
            mCurrent.y = mSearchRegion.bottom;
            advance();
            if (end) { gotoEnd();}
        }

        /* increment operators */
        circular_iterator operator++(int)
        {
            circular_iterator it = *this;
            advance();
            return it;
        }

        circular_iterator operator++()
        {
            advance();
            return *this;
        }

        /* getters */
        T& operator*()
        {
            return *(mHash->mHashMap.at(mCurrent));
        }

        T* operator&()
        {
            return mHash->mHashMap.at(mCurrent);
        }

        Point location()
        {
            return mCurrent;
        }

        /* equality of iterators - must compare compatible iterators */
        bool operator!=(const circular_iterator& other) const
        {
            if (mCenter != other.mCenter || mRadius != other.mRadius)
            {
                throw std::invalid_argument("comparison between"
                    " incompatible circular iterators");
            }

            return mCurrent != other.mCurrent;
        }

        /* goto to end of search region */
        void gotoEnd()
        {
            mCurrent.x = mSearchRegion.right + 1.0;
            mCurrent.y = mSearchRegion.top;
        }
    };

    /* iterator to full list of values */
    typedef typename std::vector<T*>::iterator iterator;

    /* constructor given minimum radius of object and bucket tolerance */
    SpatialHash(double, double = 0.01);

    /* insert & delete objects from map */
    void insert(Point, T*);
    void erase(Point, T*);

    /* update location of object */
    void update(Point, Point);

    /* get random object from map */
    T* randomValue();

    /* number of objects in map */
    int size();

    /** iterators **/

    circular_iterator begin(Point& center, double radius)
    {
        return circular_iterator(this, center, radius);
    }

    circular_iterator end(Point& center, double radius)
    {
        return circular_iterator(this, center, radius, true);
    }

    iterator begin()
    {
        return mValueList.begin();
    }

    iterator end()
    {
        return mValueList.end();
    }
};

/* constructor */
template <class T>
SpatialHash<T>::SpatialHash(double minRadius, double tol)
{
    mBucketSize = pow(2, 0.5) * minRadius / 2 - tol;
}

/* Hash function for 2D point - O(1) */
template <class T>
Point SpatialHash<T>::hash(Point pt)
{
    /* find nearest grid point */
    Point center;
    center.x = ceil((fabs(pt.x) - mBucketSize) / mBucketSize);
    center.y = ceil((fabs(pt.y) - mBucketSize) / mBucketSize);

    /* adjust if negative */    
    if (pt.x < 0) { center.x *= -1;}
    if (pt.y < 0) { center.y *= -1;}

    /* return hashed point */
    return center;
}

/* Adds key to hash map - Average Case: O(1) */
template <class T>
void SpatialHash<T>::addKey(Point pt, T* val)
{
    /* try to insert object at point, throw error on failure */
    if (!mHashMap.insert(std::pair<Point, T*>(Hash(pt), val)).second)
    {
        throw std::invalid_argument("can't add: key already mapped\n");
    }
}

/* Removes key from hash map - Average Case: O(1) */
template <class T>
void SpatialHash<T>::removeKey(Point pt)
{
    /* try to delete object from map, throw error on failure */
    if (mHashMap.erase(Hash(pt)) == 0)
    {
        throw std::invalid_argument("can't remove: key is not mapped\n");
    }
}

/*
   Inserts value into hash map - does not catch the following
   error: value already exists in mValueList but not mHashMap.
   Average Case: O(1)
 */
template <class T>
void SpatialHash<T>::insert(Point pt, T* val)
{
    mIndexMap.insert(std::pair<Point, unsigned>(Hash(pt), size()));
    mValueList.push_back(val);
    AddKey(pt, val);
}

/* Permanently deletes value - Average Case: O(1) */
template <class T>
void SpatialHash<T>::erase(Point pt, T* val)
{
    /* index in list (throws an exception if not found) */
    unsigned index = mIndexMap.at(Hash(pt));

    /* delete object for all 3 data structures */
    mValueList[index] = mValueList.back();
    mValueList.pop_back();
    mIndexMap.erase(pt);
    RemoveKey(pt);

    /* update key of moved object */
    mIndexMap.erase(mValueList[index]);
    mIndexMap.insert(std::pair<Point, unsigned>(mValueList[index], index));
}

/* Re-hash key in case it moved to new bucket - Average Case: O(1) */
template <class T>
void SpatialHash<T>::update(Point oldPt, Point newPt)
{
    /* get value */
    T* val = mHashMap.at(Hash(oldPt));

    /* remove old point, insert value at new point */
    RemoveKey(oldPt);
    AddKey(newPt, val);
}

/* Return the number of (key,value) pairs in the SpatialHash - O(1) */
template <class T>
int SpatialHash<T>::size()
{
    /* check if sizes match */
    if (mHashMap.size() != mValueList.size())
    {
        throw std::logic_error("hash map sizes out of sync\n");
    }

    /* return size */
    return mHashMap.size();
}

/* Return a random value in the SpatialHash - O(1) */
template <class T>
T* SpatialHash<T>::randomValue()
{
    /* return value at random index */
    return mValueList[floor(R::runif(0, size()))];
}

#endif

