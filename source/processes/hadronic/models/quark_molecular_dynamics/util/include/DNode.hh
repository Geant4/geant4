#ifndef _DNode_H
#define _DNode_H

#include "Definitions.hh"
#include "Boolean.h"

class DNode {
   public:
      void* object;
      DNode* prev;
      DNode* next;
      Boolean Marked;
      DNode(void* o, DNode* p = 0, DNode* n = 0) : object(o),prev(p),
                     next(n),Marked(Boolean::False) {}
};

#endif // _DNode_H

