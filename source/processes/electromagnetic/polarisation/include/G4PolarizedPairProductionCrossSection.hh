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
// $Id: G4PolarizedPairProductionCrossSection.hh,v 1.2 2006/11/17 14:14:19 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedBremsstrahlungCrossSection
//
// Author:        Karim Laihem
//
// Creation date: 15.05.2005
// 20-08-06  updated interface to geant4.8.1 (A.Schaelicke)
// Modifications:
//
// Class Description:
//   determine the  polarization of the final state 
//   in a Bremsstrahlung scattering process employing the differential 
//   cross section by Olsen and Maximon
//
#ifndef G4PolarizedPairProductionCrossSection_h
#define G4PolarizedPairProductionCrossSection_h 1

#include "G4StokesVector.hh"
#include "G4VPolarizedCrossSection.hh"
#include "G4RotationMatrix.hh"

class G4PolarizedGammaConversionModel;

class G4PolarizedPairProductionCrossSection : public G4VPolarizedCrossSection
{
 public:
  G4PolarizedPairProductionCrossSection(G4PolarizedGammaConversionModel * model);
  virtual void Initialize(G4double eps, G4double X, G4double phi,
			  const G4StokesVector & p0,
			  const G4StokesVector & p1,
			  G4int flag=0); 
  virtual G4double XSection(const G4StokesVector & pol2,
			    const G4StokesVector & pol3); 

  // return expected mean polarisation
  G4StokesVector GetPol2();  // electron/positron
  G4StokesVector GetPol3();  // photon
 private:
  
  G4PolarizedGammaConversionModel * theModel;

  G4StokesVector  theFinalElectronPolarization;
  G4StokesVector  theFinalPositronPolarization;

  void InitializeMe();

  static G4bool scrnInitialized;
  static G4double SCRN [3][20];  // screening function lookup table;
};

#endif
