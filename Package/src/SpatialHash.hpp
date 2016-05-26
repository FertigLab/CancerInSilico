#ifndef SPATIAL_HASH_HPP
#define SPATIAL_HASH

#include "Cell.hpp"

#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>

typedef struct point {

  double x;
  double y;

  bool operator<(const struct point& other) const {
    return pow(x,2) + pow(y,2) < pow(other.x,2) + pow(other.y,2);
  }

} Point;

class BucketIterator;

class SpatialHash {

private:

  std::map<Point, std::vector<Cell*> > m_buckets;
  double m_bucket_size;

  Point GetBucket(Cell*);
  Point GetBucket(Point);

public:

  friend class BucketIterator;

  SpatialHash(std::vector<Cell*>&);
  ~SpatialHash();

  void PlaceInBucket(Cell*);
  Cell* GetRandomCell();
  int size();

  BucketIterator getCircularIterator(Cell*,double);

};

class BucketIterator {

private:

  SpatialHash* hash_map;
  Point bucket;
  int index;
  std::vector<Point> ids;

public:

  BucketIterator(SpatialHash* hash, Point pt, double radius) {

    hash_map = hash;
    double y, x = pt.x - radius;

    while (x < pt.x + radius + hash->m_bucket_size) {

      y = pt.y - radius;

      while (y < pt.y + radius + hash->m_bucket_size) {

        Point p = {x, y};
        p = hash_map->GetBucket(p);
        if (hash_map->m_buckets.count(p) > 0) {
          ids.push_back(p);
        }
        y += hash->m_bucket_size;

      }

      x += hash->m_bucket_size;

    }
          
    if (!ids.empty()) {bucket = ids.back();}
    index = -1;

  }

  bool Next() {

    if (ids.empty()) {

      return false;
  
    } else if (index < hash_map->m_buckets[bucket].size() - 1) {

      ++index;
      return true;

    } else if (ids.size() == 1) {

      return false;      

    } else {

      do {

        index = 0;
        ids.pop_back();

        if (ids.empty()) {

          return false;

        }

        bucket = ids.back();

      } while (hash_map->m_buckets[bucket].empty());

      return true;

    }

  }

  Cell* getCell() {

    return hash_map->m_buckets[bucket][index];
    
  }

};

#endif
