#include <iostream>
#include "../src/G4BertiniCascade.cc" 

// test method for G4BertiniCascade::big7()

main() {

  G4bool analyze = true; // set false for verbose mode

  if (!analyze) cout << "Testing G4BertiniCascade::big7(G4double, G4int)" << endl;

  G4BertiniCascade b1;
  G4BertiniCascade b2;
  G4BertiniCascade b3;
  G4int ii=5;
  G4double d=0.3;
  for (G4int i=1; i<10000; i++){
    cout << b1.big7(d,ii) <<endl;
  }

 if (!analyze) cout << "Testing done" << endl;
 
  return 0;
}


