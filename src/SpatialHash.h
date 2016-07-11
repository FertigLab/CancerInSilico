// [[Rcpp::depends(BH)]]

#ifndef SPATIAL_HASH_HPP
#define SPATIAL_HASH_HPP

#include <R.h>
#include <boost/unordered_map.hpp>
#include <cmath>
#include <vector>

#include "Point.h"
#include "ExceptionHandling.h"

template <class T>
class SpatialHash {

friend class TestSpatialHash;
friend class TestCellPopulation;

public:

    class circular_iterator {

    public:

        circular_iterator(SpatialHash<T>* hash, Point& center, double radius, bool end = false) {
    
			m_hash = hash;
			m_center = center;
			m_radius = radius;
            construct_box(radius + 4.0);

            m_current.x = m_box.left - m_hash->m_bucket_size;
            m_current.y = m_box.bottom;
			advance_to_next();
			if (end) { goto_end();}

        }

        circular_iterator operator++(int junk) {

            circular_iterator it = *this;
            advance_to_next();
            return it;

        }

        circular_iterator operator++() {

            advance_to_next();
            return *this;

        }

        T& operator*() {
            
            return *(m_hash->m_hash_map[m_current]);

        }

		T* operator&() {

			return m_hash->m_hash_map[m_current];

		}

        bool operator==(const circular_iterator& other) const {

            return other.m_current == m_current;

        }

        bool operator!=(const circular_iterator& other) const {

			if (!(m_center == other.m_center && m_radius == other.m_radius)) {

				RCPP_STOP_TRACE("comparison between incompatible circular iterators");

			} 
	        
			return other.m_current != m_current;

        }

        void goto_end() {
        
            m_current.x = m_box.right + m_hash->m_bucket_size;
            m_current.y = m_box.top;

        }

		Point& location() {

			return m_current;
	
		}
    
    private:
    
        struct box {

            double left, right, top, bottom;

        } m_box;

        SpatialHash<T>* m_hash;
        Point m_current;
		Point m_center;
		double m_radius;

        void construct_box(double radius) {
            
            m_center.x -= radius;
            m_box.left = m_hash->Hash(m_center).x;

            m_center.x += 2 * radius;
            m_box.right = m_hash->Hash(m_center).x;            
            
            m_center.y -= radius;
            m_box.bottom = m_hash->Hash(m_center).y;            

            m_center.y += 2 * radius;
            m_box.top = m_hash->Hash(m_center).y;            
			
			m_center.x -= radius;
			m_center.y -= radius;

        }

        void advance_to_next() {

			do {
            	
				if (m_current.x > m_box.right + m_hash->m_bucket_size/2) {

					break;
				}

				m_current.y -= m_hash->m_bucket_size;
            	if (m_current.y < m_box.bottom) {
                
                	m_current.y = m_box.top;
                	m_current.x += m_hash->m_bucket_size;

            	}
				m_current = m_hash->Hash(m_current);
	
			} while (m_hash->m_hash_map.count(m_current) == 0
						|| m_current == m_hash->Hash(m_center));

        }

    };

    class full_iterator {

    public:

		full_iterator(typename std::vector<T*>::iterator iter) {

            m_iter = iter;

        }

        full_iterator operator++(int junk) {

            full_iterator ret_iter = *this;
            ++m_iter;
			
            return ret_iter;

        }

        full_iterator operator++() {

            ++m_iter;
            return *this;

        }

        T& operator*() {
            
            return *(*m_iter);

        }

		T* operator&() {

			return *m_iter;

		}

        bool operator==(const full_iterator& other) {

            return other.m_iter == m_iter;

        }

        bool operator!=(const full_iterator& other) {

            return other.m_iter != m_iter;

        }

    private:

        typename std::vector<T*>::iterator m_iter;

    };

    SpatialHash() {}
    SpatialHash(double);

    void Insert(Point, T*);
    void Delete(Point, T*);
    void Update(Point, Point);

    T* GetRandomValue();
    int size();

    circular_iterator begin(Point center, double radius) {

		return circular_iterator(this, center, radius);

	}

    circular_iterator end(Point center, double radius) {

		return circular_iterator(this, center, radius, true);

	}

    full_iterator begin() {

		size();
		return full_iterator(m_value_list.begin());

	}

    full_iterator end() {

		size();
		return full_iterator(m_value_list.end());

	}

private:

    boost::unordered_map<Point, T*, ihash, iequal_to> m_hash_map;
    std::vector<T*> m_value_list;

    double m_bucket_size;
    static double m_bucket_tol;

    Point Hash(Point);

    void RemoveKey(Point);
    void AddKey(Point, T*);

};

template <typename T>
double SpatialHash<T>::m_bucket_tol = 0.01;

template <class T>
SpatialHash<T>::SpatialHash(double min_object_radius) {

    m_bucket_size = pow(2, 0.5) * min_object_radius / 2 - m_bucket_tol;

}

// Hash function for 2D point.
// O(1)
template <class T>
Point SpatialHash<T>::Hash(Point pt) {

    Point center = Point(0,0);

    center.x = m_bucket_size * floor(fabs(pt.x) / m_bucket_size)
               + m_bucket_size / 2;

    center.y = m_bucket_size * floor(fabs(pt.y) / m_bucket_size)
               + m_bucket_size / 2;

    if (pt.x < 0) { center.x *= -1;}
    if (pt.y < 0) { center.y *= -1;}

    return center;

}

// Adds key to hash map.
// Average Case: O(1)
template <class T>
void SpatialHash<T>::AddKey(Point pt, T* val) {

    if (!m_hash_map.insert(std::pair<Point, T*>(Hash(pt), val)).second) {

        throw std::invalid_argument("can't add: key already mapped\n");

    }

}

// Removes key from hash map.
// Average Case: O(1)
template <class T>
void SpatialHash<T>::RemoveKey(Point pt) {

    if (m_hash_map.erase(Hash(pt)) == 0) {

        throw std::invalid_argument("can't remove: key is not mapped\n");

    }

}

// Inserts value into hash map - does not catch the following
// error: value already exists in m_value_list but not m_hash_map.
// Average Case: O(1)
template <class T>
void SpatialHash<T>::Insert(Point pt, T* val) {

    m_value_list.push_back(val);

    try {

        AddKey(pt, val);

    } catch (std::invalid_argument &e) {

        RCPP_STOP_TRACE("can't insert cell: key already mapped\n");

    }

}

// Permanently deletes value.
// O(N)
template <class T>
void SpatialHash<T>::Delete(Point pt, T* val) {

    typename std::vector<T*>::iterator it =
        std::find(m_value_list.begin(), m_value_list.end(), val);

    if (it != m_value_list.end()) {

		*it = m_value_list.back();
		m_value_list.pop_back();
        RemoveKey(pt);

    } else {

        RCPP_STOP_TRACE("can't delete: cell does not exis\n");

    }

}

// Re-hash key in case it moved to new bucket.
// Average Case: O(1)
template <class T>
void SpatialHash<T>::Update(Point old_pt, Point new_pt) {

	T* val = m_hash_map.at(Hash(old_pt));
    RemoveKey(old_pt);
    AddKey(new_pt, val);

}

// Return the number of (key,value) pairs in the SpatialHash.
// O(1)
template <class T>
int SpatialHash<T>::size() {

    if (m_hash_map.size() != m_value_list.size()) {

        RCPP_STOP_TRACE("hash map sizes out of sync\n");

    }

    return m_hash_map.size();

}

// Return a random value in the SpatialHash.
// O(1)
template <class T>
T* SpatialHash<T>::GetRandomValue() {

    int ndx = floor(R::runif(0, size()));
    return m_value_list[ndx];

}

#endif

