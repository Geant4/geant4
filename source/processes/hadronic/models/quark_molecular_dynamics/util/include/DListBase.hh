#ifndef _DListBase_H
#define _DListBase_H

#include "Definitions.hh"
#include "DNode.hh"
#include "Subscript.hh"
#include "Error.hh"
#include "Boolean.h"

class DListBase {
  public:
    DListBase() : head(0),tail(0),length(0) {}
    DListBase(const DListBase& l);
    ~DListBase();

    Boolean isEmpty() const { return length==0; }
    Subscript getLength() const { return length; }

    void* appendFirst(void* o);
    void* appendLast(void* o);
    void* insertAfter(void* cur,void* o);
    void* insertBefore(void* cur,void* o);
    void* removeFirst();
    void* removeLast();
    Subscript remove(void* o);
    Subscript removeAll();
    Subscript removeMarked();

    class ErrEmpty
      : public Error
    {
      virtual void writeMessage(G4std::ostream& os) const;
    };

    class IteratorBase {
      friend DListBase;
      public:
        IteratorBase(DListBase& l) : list(&l),cur(l.head) {}
        IteratorBase(const IteratorBase& i) : list(i.list),cur(i.cur) {}

        Boolean more() const { return cur != 0; }
        void* rewind() { return (cur = list->head)->object; }
        void* advance() { return (cur = cur->next) ? cur->object : 0; }
        void* current() const {
               if( cur == 0 ) {
                 throw ErrNoMore();
               }
               return cur->object;
        }
        void Mark() { cur->Marked = Boolean::True; }
      //        Boolean isMarked() const { return cur ? cur->Marked : (Boolean)Boolean::False; }

        class ErrNoMore
          : public Error
        {
          virtual void writeMessage(G4std::ostream& os) const;
        };

      private:
        DListBase* list;
        DNode* cur;
    };
    friend IteratorBase;

  private:
    DNode* head;
    DNode* tail;
    Subscript length;
};

#endif // _DListBase_H

