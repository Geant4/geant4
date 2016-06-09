//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4ionLowEnergyIonisation physics process -----
//                by Vladimir Ivanchenko, 6 September 1999 
//                was made on the base of G4hLowEnergyIonisation class
// ************************************************************
//  6 September 1999 V.Ivanchenko create
// ************************************************************
//
// Class Description:
// This class is OBSOLETE; the same functionality is now provided by 
// G4hLowEnergyIonisation 
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// ------------------------------------------------------------
 
#ifndef G4ionLowEnergyIonisation_h
#define G4ionLowEnergyIonisation_h 1
 
#include "G4hLowEnergyIonisation.hh"

class G4ionLowEnergyIonisation : public G4hLowEnergyIonisation
{
public: // Without description
  
  G4ionLowEnergyIonisation(const G4String& processName = "ionLowEIoni"); 
  
  ~G4ionLowEnergyIonisation();
    
  //  G4double GetIonParametrisedLoss(const G4Material* material, const G4double KinEnergy, 
  //		       const G4double DeltaRayCutNow);

  //  G4double GetIonBetheBlochLoss(const G4Material* material, const G4double KinEnergy,
  //		     const G4double DeltaRayCutNow);
			     
  G4double GetLowEnergyForParametrisation(const G4Material* material);
  
  void PrintInfoDefinition();

public: // With description

  void SetIonDefinition(G4ParticleDefinition* theIonType);
  // This method define the ion type

private:
  
  // hide assignment operator 
  G4ionLowEnergyIonisation & operator=(const G4ionLowEnergyIonisation &right);
  G4ionLowEnergyIonisation(const G4ionLowEnergyIonisation&);
  
  //  private data members ...............................

  //  G4double GetConstraints(const G4DynamicParticle *aParticle,
  //                      G4Material *aMaterial);
                                       
  G4double GetIonLossWithFluct(const G4DynamicParticle *aParticle,
                            G4Material *aMaterial,
                            G4double MeanLoss) ;

  //  G4VParticleChange* AlongStepDoIt(const G4Track& track ,const G4Step& Step);

protected:
  //  protected data members ...............................
    
};

#endif
 







