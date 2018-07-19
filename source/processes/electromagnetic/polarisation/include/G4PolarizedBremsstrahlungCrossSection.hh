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
// $Id: G4PolarizedBremsstrahlungCrossSection.hh 96114 2016-03-16 18:51:33Z gcosmo $
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedBremsstrahlungCrossSection
//
// Author:        Karim Laihem and Andreas Schaelicke
//
// Creation date: 15.05.2005
//
// Modifications:
//   15.10.07     introduced a more general framework for cross sections (AS)
//
// Class Description:
//   determine the  polarization of the final state 
//   in a Bremsstrahlung scattering process employing the differential 
//   cross section by Olsen and Maximon
//
#ifndef G4PolarizedBremsstrahlungCrossSection_h
#define G4PolarizedBremsstrahlungCrossSection_h 1

#include "G4StokesVector.hh"
#include "G4VPolarizedCrossSection.hh"




class G4PolarizedBremsstrahlungCrossSection : public G4VPolarizedCrossSection
{
 public:
  G4PolarizedBremsstrahlungCrossSection();
  virtual void Initialize(G4double eps, G4double X, G4double phi,
			  const G4StokesVector & p0,
			  const G4StokesVector & p1,
			  G4int flag=0) override;
  virtual G4double XSection(const G4StokesVector & pol2,
			    const G4StokesVector & pol3) override; 

  // return expected mean polarisation
  G4StokesVector GetPol2() override;  // electron/positron
  G4StokesVector GetPol3() override;  // photon
 private:

  G4StokesVector  theFinalLeptonPolarization;
  G4StokesVector  theFinalGammaPolarization;

  void InitializeMe();

  static G4bool scrnInitialized;
  static G4double SCRN [3][20];  // screening function lookup table;
};

#endif
