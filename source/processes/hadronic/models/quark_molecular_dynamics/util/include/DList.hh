#ifndef _DList_H
#define _DList_H

#ifdef IS_GCC
#pragma interface
#endif
#include "heapwr.hh"
#include "DListBase.hh"

template<class T>
class DList
  : private DListBase
{
  public:
    DListBase::isEmpty;
    DListBase::getLength;

    T* appendFirst(T* o) { return (T*) DListBase::appendFirst(o); }
    T* appendLast(T* o) { return (T*) DListBase::appendLast(o); }
    T* insertAfter(T* cur,T* o) { return (T*) DListBase::insertAfter(cur,o); }
    T* insertBefore(T* cur,T* o) { return (T*) DListBase::insertBefore(cur,o); }
    T* removeFirst() { return (T*) DListBase::removeFirst(); }
    T* removeLast() { return (T*) DListBase::removeLast(); }
    Subscript remove(T* o) { return DListBase::remove(o); }
    Subscript removeAll() { return DListBase::removeAll(); }
    Subscript removeAndDeleteAll()
    {
      Subscript count = 0;
      while( !isEmpty() ) {
        T* x = removeLast();
	if (x)
	  delete x;
        count++;
      }
      return count;
    }

    Subscript removeMarked() { return DListBase::removeMarked(); }

    class Iterator
      : public DListBase::IteratorBase
    {
      public:

        Iterator(DList<T>& l) : DListBase::IteratorBase(l) {}
        T* current() const {
           return (T*) DListBase::IteratorBase::current();
        }
        operator T*() {
           return (T*) DListBase::IteratorBase::current();
        }
        T* operator->() const {
           return (T*) DListBase::IteratorBase::current();
        }
    };
    friend Iterator;

};

#endif // _DList_H

