#include <iostream>
#include "../src/G4BertiniCascade.cc" 

// test mothod GBertiniCascade::mud()

main() {

  bool analyze = true; // set false for verbose mode

  if (!analyze) cout << "Testing G4BertiniCascade::mud()" << endl;

  G4BertiniCascade b1;
  G4BertiniCascade b2;
  G4BertiniCascade b3;

  for (int i=1; i<10000; i++){
    cout << b1.mud() << " " << b2.mud() << " " << b3.mud() <<endl;
  }

 if (!analyze) cout << "Testing done" << endl;
 
  return 0;
}


