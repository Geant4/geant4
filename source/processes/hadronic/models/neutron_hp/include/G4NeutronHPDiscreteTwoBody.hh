// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPDiscreteTwoBody.hh,v 1.2 1999-06-29 18:43:53 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPDiscreteTwoBody_h
#define G4NeutronHPDiscreteTwoBody_h 1

#include "G4ios.hh"
#include <fstream.h>
#include "globals.hh"
#include "G4VNeutronHPEnergyAngular.hh"
#include "G4NeutronHPLegendreTable.hh"
#include "G4NeutronHPInterpolator.hh"
#include "G4InterpolationManager.hh"

class G4NeutronHPDiscreteTwoBody : public G4VNeutronHPEnergyAngular
{
  public:
  
  G4NeutronHPDiscreteTwoBody();
  ~G4NeutronHPDiscreteTwoBody();
  
  void Init(ifstream & aDataFile)
  {
    aDataFile >> nEnergy;
    theManager.Init(aDataFile);
    theCoeff = new G4NeutronHPLegendreTable[nEnergy];
    for(G4int i=0; i<nEnergy; i++)
    {
      G4double energy;
      G4int aRep, nCoeff;
      aDataFile >> energy >> aRep >> nCoeff;
      energy*=eV;
      G4int nPoints=nCoeff;
      if(aRep>0) nPoints*=2;
      theCoeff[i].Init(energy, nPoints);
      theCoeff[i].SetRepresentation(aRep);
      for(G4int ii=0; ii<nPoints; ii++)
      {
        G4double y;
        aDataFile >> y;
        theCoeff[i].SetCoeff(ii, y);
      }
    }
  }
    
  G4ReactionProduct * Sample(G4double anEnergy, G4double massCode, G4double mass);
  G4double MeanEnergyOfThisInteraction() { return -1; }
  
  private:
  
  G4int nEnergy;
  G4InterpolationManager theManager; // knows the interpolation between stores
  G4NeutronHPLegendreTable * theCoeff;
    
  private:
  
  G4NeutronHPInterpolator theInt;

};
#endif
