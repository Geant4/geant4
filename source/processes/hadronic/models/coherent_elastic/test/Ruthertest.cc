// J.P. Wellisch, X-mas 2002.
#include "G4Rutherford.hh"
#include "G4ParticleTable.hh"

int main()
{
  G4Rutherford theR;
  G4ParticleDefinition * theP 
     = G4ParticleTable::GetParticleTable()->FindIon(4,8,0,4);
  G4Nucleus theN(28, 14);
  int i=0;
  while (i++<1000000) 
  {
    if(i==1000*(i/1000)) G4cerr << "Event # "<<i<<G4endl;
    G4cout << theR.Apply(theP, theN)<<G4endl;
  }
}
