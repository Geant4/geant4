// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4ionLowEnergyIonisation physics process -----
//                by Vladimir Ivanchenko, 6 September 1999 
//                was made on the base of G4hLowEnergyIonisation class
// ************************************************************
// Class Description:
// It is the extention of the ionisation process for the slow 
// charged ions. Physics model is described in CERN-OPEN-99-XXX.
// Class Description - End
// ************************************************************
//  6 September 1999 V.Ivanchenko create
// ------------------------------------------------------------
 
#ifndef G4ionLowEnergyIonisation_h
#define G4ionLowEnergyIonisation_h 1
 
#include "G4hLowEnergyIonisation.hh"

class G4ionLowEnergyIonisation : public G4hLowEnergyIonisation
{
public: // Without description
  
  G4ionLowEnergyIonisation(const G4String& processName = "ionLowEIoni"); 
  
  ~G4ionLowEnergyIonisation();
    
  G4double GetIonParametrisedLoss(const G4Material* material, const G4double KinEnergy, 
			       const G4double DeltaRayCutNow);

  G4double GetIonBetheBlochLoss(const G4Material* material, const G4double KinEnergy,
			     const G4double DeltaRayCutNow);

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

  G4double GetConstraints(const G4DynamicParticle *aParticle,
                          G4Material *aMaterial);
                                       
  G4double GetIonLossWithFluct(const G4DynamicParticle *aParticle,
                            G4Material *aMaterial,
                            G4double MeanLoss) ;

  G4VParticleChange* AlongStepDoIt(const G4Track& track ,const G4Step& Step);

  G4ParticleDefinition* theIon;

protected:
  //  protected data members ...............................


    
};

//#include "G4ionLowEnergyIonisation.icc"

#endif
 







