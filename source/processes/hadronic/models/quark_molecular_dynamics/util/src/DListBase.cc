#include <iostream.h>

#include "DListBase.hh"
#include "Boolean.h"
#include "heapwr.hh"

DListBase::DListBase(const DListBase& l)
  : head(0),
    tail(0),
    length(0)
{
  for(DNode* node = l.head; node; node = node->next) {
    appendLast(node->object);
  }
}

DListBase::~DListBase() {
  while( !isEmpty() ) {
    removeLast();
  }
}

void*
DListBase::appendFirst(void* o) {
  if( isEmpty() ) {
    tail = head = NEW DNode(o, 0, 0);
  }
  else {
    head = head->prev = NEW DNode(o, 0, head);
  }
  length++;
  return head->object;
}

void*
DListBase::appendLast(void* o) {
  if( isEmpty() ) {
    tail = head = NEW DNode(o, 0, 0);
  }
  else {
    tail = tail->next = NEW DNode(o, tail, 0);
  }
  length++;
  return tail->object;
}

void*
DListBase::insertAfter(void* obj,void* o) {
  DNode* cur = head;
  if ( isEmpty() )
     cur = tail = head = NEW DNode(o,0,0);
  else {
  while (cur != tail && cur->object != obj)
     cur = cur->next;
  if (cur->object == obj) {
     DNode* nxt = cur->next;
     cur = cur->next = NEW DNode(o,cur,nxt);
     if (nxt)
        nxt->prev = cur;
     else
        tail = cur;
  }
  else
     throw IteratorBase::ErrNoMore();
  }
  length++;
  return cur->object;
}

void*
DListBase::insertBefore(void* obj,void* o) {
  DNode* cur = head;
  if ( isEmpty() )
     cur = tail = head = NEW DNode(o,0,0);
  else {
  while (cur != tail && cur->object != obj)
     cur = cur->next;
  if (cur->object == obj) {
     DNode* prv = cur->prev;
     cur = cur->prev = NEW DNode(o,prv,cur);
     if (prv)
        prv->next = cur;
     else
        head = cur;
  }
  else
     throw IteratorBase::ErrNoMore();
  }
  length++;
  return cur->object;
}

void*
DListBase::removeFirst() {
  if( isEmpty()) {
    Throw(ErrEmpty());
  }
  void* rem = head->object;
  if(head == tail) {
    delete head;
    head = tail = 0;
  }
  else {
    DNode* node = head;
    head = head->next;
    head->prev = 0;
    delete node;
  }
  length--;
  return rem;
}

void*
DListBase::removeLast() {
  if( isEmpty()) {
    Throw(ErrEmpty());
  }
  void* rem = tail->object;
  if(head == tail) {
    delete head;
    head = tail = 0;
  }
  else {
    DNode* node = tail;
    tail = tail->prev;
    tail->next = 0;
    delete node;
  }
  length--;
  return rem;
}

Subscript
DListBase::remove(void* o) {
  int count = 0;
  DNode* node = head;
  while( node ) {
    if( node->object == o ) {
      if( node == head ) {
	removeFirst();
	node = head;
      }
      else if( node == tail ) {
	removeLast();
	node = 0;
      }
      else {
	DNode* rem = node;
	(rem->prev)->next = rem->next;
	(rem->next)->prev = rem->prev;
	node = rem->next;
	delete rem;
	length--;
      }
      count++;
    }
    else
      node = node->next;
  }
  return count;
}

Subscript
DListBase::removeAll() {
  int count = 0;
  while( !isEmpty() ) {
    ++count;
    removeLast();
  }
  return count;
}

Subscript
DListBase::removeMarked() {
  int count = 0;
  DNode* node = head;
  while( node ) {
    if( node->Marked ) {
      if( node == head ) {
	removeFirst();
	node = head;
      }
      else if( node == tail ) {
	removeLast();
	node = 0;
      }
      else {
	DNode* rem = node;
	(rem->prev)->next = rem->next;
	(rem->next)->prev = rem->prev;
	node = rem->next;
	delete rem;
	length--;
      }
      count++;
    }
    else
      node = node->next;
  }
  return count;
}

void
DListBase::ErrEmpty::writeMessage(ostream& os) const {
  os << "in DListBase::remove()" << endl;
  os << "remove from empty list";
}

void
DListBase::IteratorBase::ErrNoMore::writeMessage(ostream& os) const {
  os << "in DListBase::IteratorBase::current()" << endl;
  os << "no more items in list";
}


