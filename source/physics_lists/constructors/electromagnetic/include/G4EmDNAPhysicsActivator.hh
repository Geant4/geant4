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
// $Id: G4EmDNAPhysicsActivator.hh 66704 2013-01-10 18:20:17Z gunter $

#ifndef G4EmDNAPhysicsActivator_h
#define G4EmDNAPhysicsActivator_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EmParameters;
class G4ProcessManager;
class G4LowECapture;

class G4EmDNAPhysicsActivator : public G4VPhysicsConstructor
{
public:

  G4EmDNAPhysicsActivator(G4int ver = 1);

  virtual ~G4EmDNAPhysicsActivator();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:

  void AddElectronModels0(const G4String& region, G4LowECapture* ecap, 
			  G4bool emsc, G4double elowest, G4double elimel);
 
  void AddProtonModels0(const G4String& region, G4bool pmsc, 
			G4double elimel, G4double pminbb, G4double pmax);

  void AddHeliumModels0(const G4String& region, G4bool a1msc, G4bool a2msc, 
			G4double elimel, G4double pminbb, G4double pmax);

  void AddGenericIonModels0(const G4String& region, G4double pminbb);

  void DeactivateNuclearStopping(G4ProcessManager*, G4double elimel);
 
  G4bool HasMsc(G4ProcessManager*) const;

  G4bool IsVerbose() const;

  G4int  verbose;
  G4EmParameters* theParameters;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif






