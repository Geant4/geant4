#include <iostream>
#include "../src/G4Cascade.cc" 

// test method for G4Cascade::gtiso()

main() {

  bool analyze = true; // set false for verbose mode

  if (!analyze) cout << "Testing G4Cascade::gtiso()" << endl;

  G4Cascade cascade;
  G4double a, b, c;
  
  for (int i=0; i<10000; i++){
    cascade.gtiso(a, b, c);
    //cout << a << "\t" << b << "\t" << c <<endl;
    cout << a << endl;
  }

 if (!analyze) cout << "Testing done" << endl;
 
  return 0;
}


