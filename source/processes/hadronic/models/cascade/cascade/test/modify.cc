#include <fstream>
#include "../src/G4Cascade.cc" 

// test method for G4double G4Cascade::modify()

//lets shoot a proton with energy 0.1 MeV to a perpendicular magnetic field

int main() {

  G4Cascade cascade;
  G4double energy = 0.1, u=1, v=0, w=0;
  G4int i;

  cout <<u<<" "<<v<<endl;

  for(i=0; i < 1000; i++){
    cascade.modify(u,v,w,energy,1,0);
    cout << u <<" "<< v <<endl;
  }

  return 0;
}


