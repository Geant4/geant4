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
// $Id: G4EmProcessOptions.hh,v 1.3 2005/03/28 23:08:18 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmProcessOptions
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 27.02.2004
//
// Modifications:
//
//
// Class Description:
//
// Provide options for EM processes

// -------------------------------------------------------------------
//

#ifndef G4EmProcessOptions_h
#define G4EmProcessOptions_h 1

#include <vector>
#include "globals.hh"

class G4LossTableManager;
class G4Region;

class G4EmProcessOptions
{

public:

  G4EmProcessOptions();

  ~G4EmProcessOptions();

  void SetLossFluctuations(G4bool val);

  void SetSubCutoff(G4bool val);

  void SetIntegral(G4bool val);

  void SetMinSubRange(G4double val);

  void SetMinEnergy(G4double val);

  void SetMaxEnergy(G4double val);

  void SetMaxEnergyForPreciseRange(G4double val);

  void SetMaxEnergyForMuons(G4double val);

  void SetDEDXBinning(G4int val);

  void SetDEDXBinningForPreciseRange(G4int val);

  void SetLambdaBinning(G4int val);

  void SetStepLimits(G4double v1, G4double v2);

  void SetRandomStep(G4bool val);

  void SetApplyCuts(G4bool val);

  void SetBuildPreciseRange(G4bool val);

  void SetVerbose(G4int val, const G4String& name = "all");

  void SetLambdaFactor(G4double val);

  void ActivateDeexcitation(G4bool val, const G4Region* r = 0);

private:

  G4EmProcessOptions & operator=(const  G4EmProcessOptions &right);
  G4EmProcessOptions(const  G4EmProcessOptions&);

  G4LossTableManager* theManager;

};

//....oooOO0OOooo.......oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
