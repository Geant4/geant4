// This code implementation is the intellectual property of
// the GEANT4 collaboration.
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
//      History: New Implememtation
//    
//      ---------- G4QAOLowEnergyLoss physics process -------
//                  by Stephane Chauvie, 21 May 2000 
//
// Class Description:
// Quantal Harmonic Oscillator Model for energy loss
// of slow antiprotons 
// Class Description - End

// ************************************************************
// ------------------------------------------------------------

 
#ifndef G4QAOLowEnergyLoss_hh
#define G4QAOLowEnergyLoss_hh 1

#include "G4VhEnergyLossModel.hh"
#include "globals.hh"

class G4Material;
class G4DynamicParticle;
class G4ParticleDefinition;

class G4QAOLowEnergyLoss : public G4VhEnergyLossModel
{
public: 
  
  G4QAOLowEnergyLoss(); 
  
  ~G4QAOLowEnergyLoss();
  
  virtual G4double LowEnergyLimit() const {return 0;};
  
  virtual G4double HighEnergyLimit() const {return 0;} ;
  
  virtual G4bool IsInCharge(G4double energy, 
			    const G4ParticleDefinition* particleDefinition,
			    const G4Material*) const;

  virtual G4double EnergyLoss(const G4DynamicParticle* particle,
			      const G4Material* material) const;
 
private:
  
  // hide assignment operator 
  G4QAOLowEnergyLoss & operator=( G4QAOLowEnergyLoss &right);
  G4QAOLowEnergyLoss( G4QAOLowEnergyLoss&);
  
private:
  
  //set shell for defined material
  void SetShellMaterial(G4String materialsymbol) ;
  
  // calculate stopping number for L's term
  G4double GetL0(G4double _enorm0) ;
  G4double GetL1(G4double _enorm1) ;
  G4double GetL2(G4double _enorm2) ;
  
private:
  
  G4int nbofshell;
  G4double* shellenergy;
  G4double* shellstrength;
  
  //  variable for calculation of stopping number of L's term
  static G4double L0[67][2];
  static G4double L1[22][2];
  static G4double L2[14][2];
  
  //  material avaliable
  static G4String materialavailable[6];
  
  // materials shell energy and oscillator strenghts
  static  G4double al_en[3];
  static  G4double al_st[3];
  static  G4double si_en[3];
  static  G4double si_st[3];
  static  G4double cu_en[4];
  static  G4double cu_st[4];
  static  G4double ta_en[6];
  static  G4double ta_st[6];
  static  G4double au_en[6];
  static  G4double au_st[6];
  static  G4double pt_en[6];
  static  G4double pt_st[6];
  
}; 

//#include "G4hLowEnergyIonisation.icc"

#endif
 







