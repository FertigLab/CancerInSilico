#include <Rcpp.h>
#include <stdexcept>
#include <cmath>

#include "SpatialHash.hpp"

SpatialHash::SpatialHash(double min_cell_radius) {

  m_bucket_tol = 0.01;  
  m_bucket_size = pow(2,0.5) * min_cell_radius / 2 - m_bucket_tol;
  
}

Point SpatialHash::Hash(Cell* cell) {

  Point pt;  
  pt.x = cell->GetCoord().first;
  pt.y = cell->GetCoord().second;
    
  return Hash(pt);
  
}

Point SpatialHash::Hash(Point a) {

  Point center;

  center.x = m_bucket_size * floor(fabs(a.x) / m_bucket_size)
              + m_bucket_size / 2;
  if (a.x < 0) {center.x *= -1;}

  center.y = m_bucket_size * floor(fabs(a.y) / m_bucket_size)
              + m_bucket_size / 2;
  if (a.y < 0) {center.y *= -1;}

  return center;

}

//adds key to hash map
void SpatialHash::AddKey(Cell* cell) {

  Point key_val = Hash(cell);

  if (m_hash_map.count(key_val) == 0) {

    m_hash_map.insert(std::make_pair<Point, Cell*>(key_val, cell));

  } else {

    throw std::invalid_argument("can't add: key already mapped");

  }

}
  
//removes key from hash map
void SpatialHash::RemoveKey(Cell* cell) {

  Point key_val = Hash(cell);
  
  if (m_hash_map.count(key_val) > 0) {

    m_hash_map.erase(key_val);
    
  } else {

    throw std::invalid_argument("can't remove: key is not mapped");

  }
    
}

/** inserts cell into hash map
-does not catch the following error: cell already 
exists in m_cell_list but not m_hash_map
*/
void SpatialHash::Insert(Cell* cell) {

  m_cell_list.push_back(cell);
  AddKey(cell);

}
  
//permanently deletes cell
void SpatialHash::Delete(Cell* cell) {

  std::vector<Cell*>::iterator it = 
    std::find(m_cell_list.begin(), m_cell_list.end(), cell);
  
  if (it != m_cell_list.end()) {

    m_cell_list.erase(it);
    RemoveKey(cell);

  } else {

    throw std::invalid_argument("can't delete: cell does not exist");

  }

}

//re-hash cell in case it moved to new region
void SpatialHash::Update(Cell& orig_cell, Cell& new_cell) {

  RemoveKey(&orig_cell);
  AddKey(&new_cell);

}

int SpatialHash::size() {
  
  if (m_hash_map.size() != m_cell_list.size()) {

    throw std::runtime_error("hash map sizes out of sync");

  }

  return m_hash_map.size();

}

Cell* SpatialHash::GetRandomCell() {

  int ndx = floor(R::runif(0,size()));
  return m_cell_list[ndx];

}

SpatialIterator SpatialHash::getCircularIterator(Cell* cell, double radius) {

  Point pt;
  pt.x = cell->GetCoord().first;
  pt.y = cell->GetCoord().second;

  return SpatialIterator(this, pt, radius);

}

SpatialIterator SpatialHash::getFullIterator() {

  return SpatialIterator(this);

}

