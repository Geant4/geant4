//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPArbitaryTab.hh"
#include "G4ios.hh"

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
             << "from in G4NeutronHPArbitaryTab::Sample"<<G4endl;
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
