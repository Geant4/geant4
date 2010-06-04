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
// $Id: G4ICRU73NoDeltaModel.hh,v 1.1 2010-06-04 10:23:31 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ICRU73NoDeltaModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 04.10.2010
//
// Modifications:
//
//
// Class Description:
//
// Implementation of G4ICRU73QOModel without delta-ray

// -------------------------------------------------------------------
//

#ifndef G4ICRU73NoDeltaModel_h
#define G4ICRU73NoDeltaModel_h 1

#include "G4ICRU73QOModel.hh"

class G4ICRU73NoDeltaModel : public G4ICRU73QOModel
{

public:

  G4ICRU73NoDeltaModel(const G4ParticleDefinition* p = 0,
		       const G4String& nam = "ICRU73QONoD");

  virtual ~G4ICRU73NoDeltaModel();

  virtual G4double ComputeDEDXPerVolume(const G4Material*,
					const G4ParticleDefinition*,
					G4double kineticEnergy,
					G4double cutEnergy);

  virtual G4double CrossSectionPerVolume(const G4Material*,
					 const G4ParticleDefinition*,
					 G4double kineticEnergy,
					 G4double cutEnergy,
					 G4double maxEnergy);
private:

  // hide assignment operator
  G4ICRU73NoDeltaModel & operator=(const  G4ICRU73NoDeltaModel &right);
  G4ICRU73NoDeltaModel(const  G4ICRU73NoDeltaModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
