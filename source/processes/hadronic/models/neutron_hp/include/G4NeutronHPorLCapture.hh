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
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

//
// Hadronic Process: High Precision low E neutron tracking
// original by H.P. Wellisch, TRIUMF, 14-Feb-97
// Builds and has the Cross-section data for one material.
 
#ifndef G4NeutronHPorLCapture_h
#define G4NeutronHPorLCapture_h 1

// Class Description
// Final state production model for a high precision (based on evaluated data
// libraries) description of neutron elastic scattering below 20 MeV; 
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

#include "G4NeutronHPorLCaptureData.hh"

#include "globals.hh"
#include "G4NeutronHPChannel.hh"
#include "G4HadronicInteraction.hh"

#include <set>

class G4NeutronHPorLCapture : public G4HadronicInteraction
{
  public: 
  
  G4NeutronHPorLCapture();
  
  ~G4NeutronHPorLCapture();
  
  G4HadFinalState * ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus);

  

  G4int GetNiso() {return theCapture[0].GetNiso();}
  private:
  
  G4double * xSec;
  G4NeutronHPChannel * theCapture;
  G4String dirName;
  G4int numEle;

   public: 
      G4bool IsThisElementOK ( G4String );
   private:
     std::set< G4String > unavailable_elements;

   public: 
      G4VCrossSectionDataSet* GiveXSectionDataSet() { return theDataSet; }; 
   private:
      G4NeutronHPorLCaptureData* theDataSet;
      void createXSectionDataSet(); 
};

#endif
