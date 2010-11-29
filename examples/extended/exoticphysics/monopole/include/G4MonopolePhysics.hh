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
// $Id: G4MonopolePhysics.hh,v 1.4 2010-11-29 15:14:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  12.07.10  S.Burdin (changed the magnetic and electric charge variables from integer to double)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4MonopolePhysics_h
#define G4MonopolePhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MonopolePhysicsMessenger;
class G4Monopole;

class G4MonopolePhysics : public G4VPhysicsConstructor
{
public:

  G4MonopolePhysics(const G4String& nam = "Monopole Physics");
  ~G4MonopolePhysics();

  // This method is dummy for physics
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();

  void SetMagneticCharge(G4double);
  void SetElectricCharge(G4double);
  void SetMonopoleMass(G4double);

private:

  // hide assignment operator
  G4MonopolePhysics & operator=(const G4MonopolePhysics &right);
  G4MonopolePhysics(const G4MonopolePhysics&);

  G4double    magCharge;
  G4double    elCharge;
  G4double monopoleMass;

  G4MonopolePhysicsMessenger*  theMessenger;
  G4Monopole* mpl;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

