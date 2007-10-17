
#ifndef TestCram_hh
#define TestCram_hh 1

#include "TObject.h"
#include "InclAblaTest.hh"

class TestCram : public InclAblaTest {

public:
  TestCram();
  ~TestCram();

  void runMe();

private:

  ClassDef(TestCram, 0) // Add documentation string here
};

#endif
