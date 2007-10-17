
#ifndef TestBarfit_hh
#define TestBarfit_hh 1

#include "TObject.h"
#include "InclAblaTest.hh"

class TestBarfit : public InclAblaTest {

public:
  TestBarfit();
  ~TestBarfit();

  void runMe();

private:

  ClassDef(TestBarfit, 0) // Add documentation string here
};

#endif
