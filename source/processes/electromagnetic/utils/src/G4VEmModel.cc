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
// $Id: G4VEmModel.cc,v 1.6 2006/06/29 19:55:17 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VEmModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 25.07.2005
//
// Modifications:
// 25.10.2005 Set default highLimit=100.TeV (V.Ivanchenko)
// 06.02.2006 add method ComputeMeanFreePath() (mma)
//
//
// Class Description:
//
// Abstract interface to energy loss models

// -------------------------------------------------------------------
//

#include "G4VEmModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VEmModel::G4VEmModel(const G4String& nam):
 lowLimit(0.1*keV), highLimit(100.0*TeV), fluc(0), name(nam), pParticleChange(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VEmModel::~G4VEmModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::CrossSectionPerVolume(const G4Material* material,
					   const G4ParticleDefinition* p,
					         G4double ekin,
					         G4double emin,
                                                 G4double emax)
{
  G4double cross = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
  size_t nelm = material->GetNumberOfElements();
  for (size_t i=0; i<nelm; i++) {
    const G4Element* elm = (*theElementVector)[i];
    cross += theAtomNumDensityVector[i]*
      ComputeCrossSectionPerAtom(p,ekin,elm->GetZ(),elm->GetN(),emin,emax);
    xsec[i] = cross;
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::ComputeMeanFreePath(const G4ParticleDefinition* p,
					       G4double ekin,
					 const G4Material* material,     
					       G4double emin,
                                               G4double emax)
{
  G4double mfp = DBL_MAX;
  G4double cross = CrossSectionPerVolume(material,p,ekin,emin,emax);
  if (cross >0.) mfp = 1./cross;
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


