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
// $Id: G4hBremsstrahlungModel.hh 97391 2016-06-02 10:08:45Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4hBremsstrahlungModel
//
// Author:        Vladimir Ivanchenko on base of G4MuBremsstrahlungModel
// 
// Creation date: 28.02.2008
//
// Modifications:
//
//

//
// Class Description:
//
// Implementation of energy loss for gamma emission by hadrons

// -------------------------------------------------------------------
//

#ifndef G4hBremsstrahlungModel_h
#define G4hBremsstrahlungModel_h 1

#include "G4MuBremsstrahlungModel.hh"

class G4hBremsstrahlungModel : public G4MuBremsstrahlungModel
{

public:

  explicit G4hBremsstrahlungModel(const G4ParticleDefinition* p = nullptr,
			 const G4String& nam = "hBrem");

  virtual ~G4hBremsstrahlungModel();

protected:

  virtual G4double ComputeDMicroscopicCrossSection(G4double tkin,
						   G4double Z,
						   G4double gammaEnergy) override;
private:

  // hide assignment operator
  G4hBremsstrahlungModel & 
    operator=(const  G4hBremsstrahlungModel &right) = delete;
  G4hBremsstrahlungModel(const  G4hBremsstrahlungModel&) = delete;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
