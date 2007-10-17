
#ifndef TestSpdef_hh
#define TestSpdef_hh 1

#include "TObject.h"
#include "InclAblaTest.hh"

class TestSpdef : public InclAblaTest {

public:
  TestSpdef();
  ~TestSpdef();

  void runMe();

private:

  ClassDef(TestSpdef, 0) // Add documentation string here
};

#endif
