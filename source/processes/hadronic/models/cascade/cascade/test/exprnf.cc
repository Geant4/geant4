#include <iostream>
#include "../src/G4Cascade.cc" 

// test method for G4Cascade::exprnf()

main() {

  bool analyze = true; // set false for verbose mode

  if (!analyze) cout << "Testing G4Cascade::exprnf()" << endl;

  G4Cascade c;
  
  for (int i=0; i<10000; i++){
    cout << c.exprnf()  <<endl;
  }

 if (!analyze) cout << "Testing done" << endl;
 
  return 0;
}


