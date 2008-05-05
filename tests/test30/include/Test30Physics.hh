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
// -------------------------------------------------------------
//      GEANT 4 class 
//
//      ---------- Test30Physics -------
//                by Vladimir Ivanchenko, 12 March 2002 
// 
//    Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Test30Physics_h
#define Test30Physics_h 1

#include "globals.hh"
#include "Test30HadronProduction.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4BinaryCascade.hh"

class G4VProcess;
class G4Material;
class G4QuasiElasticChannel;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test30Physics
{
public:

  Test30Physics();

  ~Test30Physics();

  G4VProcess* GetProcess(const G4String&, const G4String&, G4Material*);

  G4double GetNucleusMass() {return theProcess->GetMass();};

  G4ExcitationHandler* GetDeExcitation() {return theDeExcitation;};

  G4PreCompoundModel* GetPreCompound() {return thePreCompound;};

  void SetA(G4int A) {if(theProcess) theProcess->SetA(A);};

//    void setCutOnP(G4double val) {if(hkmod) hkmod->setCutOnP(val);};
//    void setCutOnPPP(G4double val) {if(hkmod) hkmod->setCutOnPPP(val);};

private:

  void Initialise();

  Test30HadronProduction* theProcess;
  G4ExcitationHandler*    theDeExcitation;
  G4PreCompoundModel*     thePreCompound;
  G4BinaryCascade*        hkmod;
  G4QuasiElasticChannel*  theQuasiElastic;
};

#endif

 


