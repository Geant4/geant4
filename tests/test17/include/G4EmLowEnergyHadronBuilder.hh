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
// $Id: G4EmLowEnergyHadronBuilder.hh,v 1.1 2004-05-26 11:39:08 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4EmLowEnergyHadronBuilder_h
#define G4EmLowEnergyHadronBuilder_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsListMessenger;
class G4hLowEnergyIonisation;

class G4EmLowEnergyHadronBuilder : public G4VPhysicsConstructor
{
public:
  G4EmLowEnergyHadronBuilder(PhysicsListMessenger* messenger,
                       const G4String& name = "EM_lowe_had");
  virtual ~G4EmLowEnergyHadronBuilder();

public:
  // This method is dummy for physics
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();
  
  // Activate Auger and Fluorescence
  void SetGammaCut(G4double);
  void SetAugerCut(G4double);
  void ActivateFluorescence(G4bool);
  void ActivateAuger(G4bool);

private:

   // hide assignment operator
  G4EmLowEnergyHadronBuilder & operator=(const G4EmLowEnergyHadronBuilder &right);
  G4EmLowEnergyHadronBuilder(const G4EmLowEnergyHadronBuilder&);

  std::vector<G4hLowEnergyIonisation*>  hio;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

