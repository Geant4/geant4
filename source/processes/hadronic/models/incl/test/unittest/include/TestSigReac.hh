
#ifndef TestSigReac_hh
#define TestSigReac_hh 1

#include "TObject.h"
#include "InclAblaTest.hh"

class TestSigReac : public InclAblaTest {

public:
  TestSigReac();
  ~TestSigReac();

  void runMe();

private:

  ClassDef(TestSigReac, 0) // Add documentation string here
};

#endif
