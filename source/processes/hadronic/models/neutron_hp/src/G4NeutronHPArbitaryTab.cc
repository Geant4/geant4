// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPArbitaryTab.hh"
#include "G4ios.hh"

  G4NeutronHPArbitaryTab::G4NeutronHPArbitaryTab()
  {
   theDistFunc = NULL;
  }
  G4NeutronHPArbitaryTab::~G4NeutronHPArbitaryTab()
  {
   if(theDistFunc!=NULL) delete [] theDistFunc;
  }
  G4double G4NeutronHPArbitaryTab::Sample(G4double anEnergy) 
  {
    G4int i;
    for(i=0;i<nDistFunc;i++)
    {
      if(anEnergy<theDistFunc[i].GetLabel()) break; // that is the energy we need
    }
    G4int low, high;
    if(i==nDistFunc) 
    {
      low = i-2;
      high = i-1;
    }
    else if(i==0)
    {
      if(nDistFunc==0)
      {
        G4cerr << "No distribution functions to sample "
             << "from in G4NeutronHPArbitaryTab::Sample"<<endl;
        G4Exception();
      } 
      else 
      {
        return theDistFunc[0].Sample();
      }
    }
    else
    {
      low = i-1;
      high = i;
    }
    theBuffer.Merge(theManager.GetScheme(low), anEnergy, 
                    theDistFunc+low, theDistFunc+high);
    return theBuffer.Sample();
  }
