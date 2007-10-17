
#ifndef TestBipol_hh
#define TestBipol_hh 1

#include "TObject.h"
#include "InclAblaTest.hh"

class TestBipol : public InclAblaTest {

public:
  TestBipol();
  ~TestBipol();

  void runMe();

private:

  ClassDef(TestBipol, 0) // Add documentation string here
};

#endif
