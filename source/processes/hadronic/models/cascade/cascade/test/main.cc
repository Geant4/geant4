#include <iostream.h>
#include "G4Cascade.cc"
 
main()
{
  G4Cascade cascade;

  cout << "Testing cascade" << endl;
  cascade.setVerboseLevel(1);
  cascade.go();

  cascade.setVerboseLevel(2);
  cascade.go();

  cout << "Testing done." << endl;
  return 0;
}
