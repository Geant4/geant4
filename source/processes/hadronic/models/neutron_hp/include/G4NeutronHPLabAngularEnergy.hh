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
//
// $Id: G4NeutronHPLabAngularEnergy.hh,v 1.5 2001-07-11 10:07:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPLabAngularEnergy_h
#define G4NeutronHPLabAngularEnergy_h 1

#include "G4ios.hh"
#include "g4std/fstream"
#include "globals.hh"
#include "G4Neutron.hh"
#include "G4NeutronHPInterpolator.hh"
#include "G4NeutronHPVector.hh"
#include "G4VNeutronHPEnergyAngular.hh"
#include "G4ReactionProduct.hh"
#include "G4InterpolationManager.hh"

class G4NeutronHPLabAngularEnergy : public G4VNeutronHPEnergyAngular
{
  public:
  
  G4NeutronHPLabAngularEnergy()
  {
    theEnergies = NULL;
    theData = NULL;
    nCosTh = NULL;
    theSecondManager = NULL;
  }
  ~G4NeutronHPLabAngularEnergy()
  {
    if(theEnergies != NULL) delete [] theEnergies;
    if(nCosTh != NULL) delete [] nCosTh;
    if(theData != NULL) 
    {
      for(G4int i=0; i<nEnergies; i++)
        delete [] theData[i];
      delete [] theData;
    }
    if(theSecondManager != NULL) delete [] theSecondManager;
  }
  
  public:
  
  void Init(G4std::ifstream & aDataFile);
     G4ReactionProduct * Sample(G4double anEnergy, G4double massCode, G4double mass);
  G4double MeanEnergyOfThisInteraction()
  {
    return currentMeanEnergy;
  }
  
  
  private:
  
  // number of incoming neutron energies
  G4int nEnergies;
  // Interpol between neutron energies
  G4InterpolationManager theManager; 
  // Incoming neutron energies
  G4double * theEnergies; 
  // number of directioncosines; parallel to theEnergies
  G4int * nCosTh; 
  // knows the interpolation between these stores
  G4InterpolationManager * theSecondManager; 
  // vectors of secondary energy, haufigkeit; parallel to theEnergies
  G4NeutronHPVector ** theData; 

  // utility interpolator
  G4NeutronHPInterpolator theInt;
  
  // cashed value of mean secondary energy in this event.
  G4double currentMeanEnergy;
};
#endif
