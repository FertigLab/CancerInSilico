#include "Grid.h"

template <class T>
Grid<T>::Grid(double width) {}

template <class T>
void Grid<T>::insert(const Point& key, const T* val) {}

template <class T>
void Grid<T>::erase(const Point& key, const T* val) {}

template <class T>
void Grid<T>::update(const Point& oldLoc, const Point& newLoc) {}

template <class T>
T* Grid<T>::randomValue() const {}

template <class T>
int Grid<T>::size() const {}

template <class T>
circular_iterator Grid<T>::begin(const Point&, double) const {}

template <class T>
circular_iterator Grid<T>::end(const Point&, double) const {}

template <class T>
iterator Grid<T>::begin(const Point&, double) const {}

template <class T>
iterator Grid<T>::end(const Point&, double) const {}

typedef typename Grid<T>::circular_iterator circIterator;

void circIterator::constructRegion(double radius) {}
void circIterator::advance(double radius) {}


