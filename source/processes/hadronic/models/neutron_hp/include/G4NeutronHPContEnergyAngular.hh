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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4NeutronHPContEnergyAngular.hh,v 1.6 2001-07-26 09:27:54 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPContEnergyAngular_h
#define G4NeutronHPContEnergyAngular_h 1

#include "G4ios.hh"
#include "g4std/fstream"
#include "globals.hh"
#include "G4VNeutronHPEnergyAngular.hh"
#include "G4NeutronHPContAngularPar.hh"
#include "G4InterpolationManager.hh"

// we will need one of these per product.

class G4NeutronHPContEnergyAngular : public G4VNeutronHPEnergyAngular
{
  public:
  
  G4NeutronHPContEnergyAngular()
  {
    theAngular = NULL;
    currentMeanEnergy = -2;
  }
  
  ~G4NeutronHPContEnergyAngular()
  {
    if(theAngular!=NULL) delete [] theAngular;
  }
  
  public:
  
  void Init(G4std::ifstream & aDataFile)
  {
    aDataFile >> theTargetCode >> theAngularRep >> theInterpolation >> nEnergy;
    theAngular = new G4NeutronHPContAngularPar[nEnergy];
    theManager.Init(aDataFile);
    for(G4int i=0; i<nEnergy; i++)
    {
      theAngular[i].Init(aDataFile);
      theAngular[i].SetInterpolation(theInterpolation);
    }
  }
G4double MeanEnergyOfThisInteraction();
G4ReactionProduct * Sample(G4double anEnergy, G4double massCode, G4double mass);
  
  private:
  
  G4double theTargetCode;
  G4int theAngularRep;
  G4int nEnergy;
  
  G4int theInterpolation;

  G4InterpolationManager theManager; // knows the interpolation between stores
  G4NeutronHPContAngularPar * theAngular;
  
  G4double currentMeanEnergy;
  
};
#endif
