// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//

#include "G4NeutronHPFieldPoint.hh"

  G4NeutronHPFieldPoint::G4NeutronHPFieldPoint()
  {
    X = 0;
    nP = 0;
    Y = NULL;
  }

G4NeutronHPFieldPoint::G4NeutronHPFieldPoint(G4int n)
  {
    nP = n;
    X = 0;
    Y = new G4double[nP];
    for (G4int i=0; i<nP; i++) Y[i]=0.;
  }
  
void G4NeutronHPFieldPoint::operator= (const G4NeutronHPFieldPoint & aSet)
  {
    if(&aSet!=this)
    {
      X = aSet.GetX();
      if(Y!=NULL) delete [] Y;
      Y = new G4double[aSet.GetDepth()];
      for(G4int i=0; i<aSet.GetDepth(); i++) Y[i] = aSet.GetY(i);
    }
  }

G4NeutronHPFieldPoint::~G4NeutronHPFieldPoint()
  {
   if(Y!=NULL) delete [] Y;
  }
    
void G4NeutronHPFieldPoint::InitY(G4int n)
  {
    nP = n;
    X=0;
    Y = new G4double[nP];
    for (G4int i=0; i<nP; i++) Y[i]=0.;
  }
