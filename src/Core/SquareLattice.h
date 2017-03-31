#ifndef CIS_SQUARE_LATTICE_H
#define CIS_SQUARE_LATTICE_H

#include <cmath>

#define ROOT_2 1.414213562

#include "Lattice.h"

template <class T>
class SquareLattice : public Lattice<T>
{
private:
    // space between grid lines
    double mGridWidth;

    // hash a point in 2D space to a square grid
    GridPoint hash(const Point<double>& pt) const
    {
        GridPoint hashedPt;
        hashedPt.x = ceil((fabs(pt.x) - mGridWidth / 2) / mGridWidth);
        hashedPt.y = ceil((fabs(pt.y) - mGridWidth / 2) / mGridWidth);

        hashedPt.x *= pt.x < 0 ? -1.0 : 1.0;
        hashedPt.y *= pt.y < 0 ? -1.0 : 1.0;

        return hashedPt;
    }

    /*********** local iterator implementation ************/
    class SquareLocalIterator : public BaseLocalIterator<T>
    {
    private:

        struct {int left, right, top, bottom;} mSearchRegion;
        GridPoint mHashedCenter;

    public:

        // needed for permission and mGridWidth
        #define SQ_LAT ((SquareLattice*) this->mLattice)

        SquareLocalIterator(Lattice<T>* lat, const Point<double>& center,
        double radius, bool end) : BaseLocalIterator<T>(lat, center, radius)
        {
            mHashedCenter = SQ_LAT->hash(center);
            int rad = floor((radius + SQ_LAT->mGridWidth * ROOT_2)
                / SQ_LAT->mGridWidth);

            mSearchRegion.left = mHashedCenter.x - rad;
            mSearchRegion.right = mHashedCenter.x + rad;
            mSearchRegion.bottom = mHashedCenter.y - rad;
            mSearchRegion.top = mHashedCenter.y + rad;

            if (end)
            {
                this->mCurrent.x = mSearchRegion.right + 1;
                this->mCurrent.y = mSearchRegion.top;
            }
            else
            {
                this->mCurrent.x = mSearchRegion.left - 1;
                this->mCurrent.y = mSearchRegion.bottom;
                this->operator++();
            }
        }

        BaseLocalIterator<T>* newCopy() const
            {return new SquareLocalIterator(*this);}

        void operator++()
        {
            do
            {
                if (this->mCurrent.x > mSearchRegion.right) {break;}
        
                this->mCurrent.y--;
                if (this->mCurrent.y < mSearchRegion.bottom)
                {
                    this->mCurrent.y = mSearchRegion.top;
                    this->mCurrent.x++;
                }
            } while (!SQ_LAT->mGrid.count(this->mCurrent));
        }
    };
    /******************************************************/

public:
    
    // constructor
    SquareLattice(double width) {mGridWidth = width;}

    // convenient typedef for notation
    typedef typename Lattice<T>::local_iterator LocIter;

    // local iterators
    virtual LocIter lbegin(const Point<double>& c, double r)
        {return LocIter(new SquareLocalIterator(this, c, r, false));}
    virtual LocIter lend(const Point<double>& c, double r)
        {return LocIter(new SquareLocalIterator(this, c, r, true));}
};

#endif
