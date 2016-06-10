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
// $Id: G4hBremsstrahlungModel.cc 74020 2013-09-19 13:38:38Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
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
// Class Description:
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4hBremsstrahlungModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"

using namespace std;

G4hBremsstrahlungModel::G4hBremsstrahlungModel(const G4ParticleDefinition* p,
					       const G4String& nam)
  : G4MuBremsstrahlungModel(p, nam)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4hBremsstrahlungModel::~G4hBremsstrahlungModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4hBremsstrahlungModel::ComputeDMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double gammaEnergy)
//  differential cross section
{
  G4double dxsection = 0.;

  if( gammaEnergy > tkin) return dxsection ;
  //  G4cout << "G4hBremsstrahlungModel m= " << mass 
  //	 << "  " << particle->GetParticleName() << G4endl;
  G4double E = tkin + mass ;
  G4double v = gammaEnergy/E ;
  G4double delta = 0.5*mass*mass*v/(E-gammaEnergy) ;
  G4double rab0=delta*sqrte ;

  G4int iz = G4int(Z);
  if(iz < 1) { iz = 1; }

  G4double z13 = 1.0/nist->GetZ13(iz);
  G4double dn  = mass*nist->GetA27(iz)/(70.*MeV);

  G4double    b = btf;
  if(1 == iz) b = bh;

  // nucleus contribution logarithm
  G4double rab1=b*z13;
  G4double fn=G4Log(rab1/(dn*(electron_mass_c2+rab0*rab1))*
              (mass+delta*(dn*sqrte-2.))) ;
  if(fn <0.) fn = 0. ;

  G4double x = 1.0 - v;
  if(particle->GetPDGSpin() != 0) { x += 0.75*v*v; }

  dxsection = coeff*x*Z*Z*fn/gammaEnergy;

  return dxsection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
