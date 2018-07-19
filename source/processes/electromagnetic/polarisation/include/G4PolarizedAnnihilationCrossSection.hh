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
// -------------------------------------------------------------------
// $Id: G4PolarizedAnnihilationCrossSection.hh 96114 2016-03-16 18:51:33Z gcosmo $
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     PolarizedAnnihilationCrossSection
//
// Author:        Andreas Schaelicke and Pavel Starovoitov
//
// Creation date: 22.03.2006
//
// Modifications:
//   15.10.07     introduced a more general framework for cross sections
//
// Class Description:
//   * calculates the differential cross section (ME squared, 
//     without phase space) incoming positron (along positive z direction) 
//     annihilations with an electron at rest 
//   * phi denotes the angle between the scattering plane and 
//     X axis of incoming partice reference frame (PRF) 
//
#ifndef G4PolarizedAnnihilationCrossSection_h
#define G4PolarizedAnnihilationCrossSection_h 1

#include "G4StokesVector.hh"
#include "G4VPolarizedCrossSection.hh"

class G4PolarizedAnnihilationCrossSection : public G4VPolarizedCrossSection
{
public:
  G4PolarizedAnnihilationCrossSection();
  virtual ~G4PolarizedAnnihilationCrossSection();
public:
  virtual void Initialize(G4double eps, G4double gamma, G4double phi, 
		  const G4StokesVector & p0,const G4StokesVector & p1,
		  G4int flag=0) override; 

  G4double DiceEpsilon(); 
  virtual G4double XSection(const G4StokesVector & pol2,
			    const G4StokesVector & pol3) override; 
  virtual G4double TotalXSection(G4double xmin, G4double xmax, 
				 G4double y,
				 const G4StokesVector & pol0,
				 const G4StokesVector & pol1) override; 

  // return expected mean polarisation
  G4StokesVector GetPol2() override;
  G4StokesVector GetPol3() override;

  virtual G4double GetXmin(G4double y) override; // minimal energy fraction in TotalXSection
  virtual G4double GetXmax(G4double y) override; // maximal energy fraction in TotalXSection

  G4double getVar(G4int );
  // test routine
  void getCoeff();
private:
  void TotalXS();
  void DefineCoefficients(const G4StokesVector & pol0,
    		          const G4StokesVector & pol1);


  G4double   polxx, polyy, polzz, polxz, polzx, polxy, polyx, polyz, polzy;

  G4double re2, diffXSFactor, totalXSFactor;
  // - unpolarised + part depending on the polarization of the intial pair
  G4double phi0;
  // - part depending on the polarization of the final positron
  G4ThreeVector phi2; 
  // - part depending on the polarization of the final electron
  G4ThreeVector phi3;
  G4double dice;
  G4double polXS, unpXS;
  G4double ISPxx, ISPyy, ISPzz, ISPnd;
};
#endif
