#include <iostream>
#include "../src/G4Cascade.cc"

main(){ 
  //test method for void G4Cascade::readh(G4double geosig[])

  bool analyze = false;

  if(!analyze){
    cout << "Testing void G4Cascade::readh(G4double geosig[])" << endl; 
    cout << "The numbers in vector geosig after the call of readh()
             (indexes 0-240):" << endl;
  }

  G4double geosig[241];
  G4Cascade c;

  for(int i=0; i<=240; i++)
    geosig[i]=0;

  c.readh(geosig); 

  for(int i=0; i <= 240; i++)
    cout << geosig[i] << endl; 
 
  if(!analyze)
    cout << "Testing done." << endl;

  return 0;
 
}
