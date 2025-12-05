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

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4Event;
class G4ParticleGun;

class G4IAEAphspReader;
class PrimaryGeneratorMessenger;
class DetectorConstruction;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  // The constructor defines a ParticleGun object, which allows
  // shooting a beam of particles through the experimental set-up.
  // It also needs a pointer to G4IAEAphspReader object in case we need to
  // read particles from an IAEAphsp file
  PrimaryGeneratorAction(const G4int threads);

  //The destructor. It deletes the ParticleGun.
  virtual ~PrimaryGeneratorAction() override;
 
  //Generates the primary event via the ParticleGun method,
  // and from the IAEA phase-space file.
  void GeneratePrimaries(G4Event* anEvent) override;

  //Get/Set methods
  inline void SetKinE(const G4double val)  { fKinE = val; };
  inline void SetDE(const G4double val)    { fDE = val; };
  inline void SetX0(const G4double val)  { fX0 = val;};
  inline void SetY0(const G4double val)  { fY0 = val;};
  inline void SetZ0(const G4double val)  { fZ0 = val;};
  inline void SetDX(const G4double val)  { fDX = val;};
  inline void SetDY(const G4double val)  { fDY = val;};
  inline void SetDZ(const G4double val)  { fDZ = val;};
  inline void SetVerbose(const G4int val) { fVerbose = val;};

  void SetIAEAphspReader(const G4String filename);

  inline G4int GetVerbose() const { return fVerbose; }
  inline G4IAEAphspReader* GetIAEAphspReader() const  {return fIAEAphspReader;}


private:

  // Phase space reader
  G4IAEAphspReader* fIAEAphspReader = nullptr;
  G4int             fThreads;
  G4String          fIAEAphspReaderName;

  G4int fVerbose;
  PrimaryGeneratorMessenger* fMessenger;
  G4ParticleGun* fParticleGun;

  G4int fCounter;
  G4double fKinE, fDE;
  G4double fX0, fY0, fZ0;
  G4double fDX, fDY, fDZ;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
