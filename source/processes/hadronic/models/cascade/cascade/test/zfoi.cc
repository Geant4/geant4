#include <fstream>
#include "../src/G4Cascade.cc" 

//test for G4Cascade::zfoi
main() {

  bool analyze = true; // set false for verbose mode
  G4int i;
  G4double param;
  G4Cascade cascade;

  if (!analyze) cout << "Testing G4Cascade::zfoi()" << endl;
    
  
  for(i=0; i < 410; i++){
    param = -20.0 + 0.1*i; 
    cout << param <<" "<< static_cast<G4int>(param)<<" "<< cascade.zfoi(param)
 << endl;
  }

  if (!analyze) cout << "Testing done" << endl;
 
  return 0;
}


