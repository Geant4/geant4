#include "g4std/iostream"
#include "G4Cascade.cc"
 
main()
{
  G4Cascade cascade;

  cout << "Testing cascade" << G4endl;
  cascade.setVerboseLevel(1);
  cascade.go();

  cascade.setVerboseLevel(2);
  cascade.go();

  cout << "Testing done." << G4endl;
  return 0;
}
