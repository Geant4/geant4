
#ifndef InclAblaF_hh
#define InclAblaF_hh 1

#include "TObject.h"

#include "functionwrapper.hh"

class InclAblaF : public TObject {

public:
  // Empty constructor
  InclAblaF();
  ~InclAblaF();

  // Static methods
  static void ribmF(float *rndm, int *ial);
  static void testrnF(double *rnumber);

  ClassDef(InclAblaF, 1)
};

#endif
