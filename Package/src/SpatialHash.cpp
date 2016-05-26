#include "SpatialHash.hpp"
#include <iostream>

SpatialHash::SpatialHash(std::vector<Cell*>& cell_list) {

  m_bucket_size = 0.5;

  std::vector<Cell*>::iterator iter = cell_list.begin();

  for (; iter != cell_list.end(); ++iter) {
    
    PlaceInBucket(*iter);

  }

}

SpatialHash::~SpatialHash() {}

Point SpatialHash::GetBucket(Cell* cell) {

  Point pt;  
  pt.x = cell->GetCoord().first;
  pt.y = cell->GetCoord().second;
    
  return GetBucket(pt);
  
}

Point SpatialHash::GetBucket(Point a) {

  Point center;

  center.x = floor(abs(a.x) / m_bucket_size) - m_bucket_size / 2;
  if (a.x < 0) {center.x *= -1;}

  center.y = floor(abs(a.y) / m_bucket_size) - m_bucket_size / 2;
  if (a.y < 0) {center.y *= -1;}

  return center;

}

void SpatialHash::PlaceInBucket(Cell* cell) {
  
  Point bucket_id = GetBucket(cell);

  if (m_buckets.count(bucket_id) == 0) {

    m_buckets[bucket_id] = std::vector<Cell*>();

  }

  m_buckets[bucket_id].push_back(cell);
  
}

int SpatialHash::size() {
  
  int sum = 0;
  std::map<Point,std::vector<Cell*> >::iterator iter = m_buckets.begin();

  for (; iter != m_buckets.end(); ++iter) {

    sum += iter->second.size();

  }

  return sum;

}

Cell* SpatialHash::GetRandomCell() {

  int ndx = rand() % size();

  std::map<Point,std::vector<Cell*> >::iterator iter = m_buckets.begin();

  for (; iter != m_buckets.end(); ++iter) {

    if (ndx < iter->second.size()) {

      return iter->second[ndx];

    } else {

      ndx -= iter->second.size();

    }

  }

}

BucketIterator SpatialHash::getCircularIterator(Cell* cell, double radius) {

  Point pt;
  pt.x = cell->GetCoord().first;
  pt.y = cell->GetCoord().second;

  return BucketIterator(this, pt, radius);

}


