#include <iostream>
#include "../src/G4Cascade.cc" 

// test method for G4double G4Cascade::xlamb(G4double x, G4double y, G4double z)

int main() {

  G4bool analyze = true; // set false for verbose mode
  G4Cascade cascade;
  G4double xtable[10000], ytable[10000], ztable[10000], result[10000];

  if (!analyze) cout << "Testing G4Cascade::xlamb()" << endl;

  
  for (G4int i=0; i < 10000; i++){
    xtable[i]=-5+0.001*i;
    ytable[i]=-5+0.001*i;
    ztable[i]=-5+0.001*i;

    result[i]=cascade.xlamb(xtable[i], ytable[i], ztable[i]);

    if(result[i]<0){
      if(!analyze)
	cout << "Test failed" << endl;
      break; 
    }
    
    cout << xtable[i] <<" "<< ytable[i] << " " << ztable[i] << " " <<
 result[i] << endl;
  }

 if (!analyze) cout << "Testing done" << endl;
 
  return 0;
}


