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
// $Id: G4PolarizedComptonCrossSection.hh 96114 2016-03-16 18:51:33Z gcosmo $
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedComptonCrossSection
//
// Author:        Andreas Schaelicke
//
// Creation date: 15.05.2005
//
// Modifications:
// 18-07-06 use newly calculated cross sections (P. Starovoitov)
// 21-08-06 update interface to geant4.8.1 (A. Schaelicke)
// 15-10-07 introduced a more general framework for cross sections (AS)
//
//
// Class Description:
//   determine the  polarization of the final state 
//   in a Compton scattering process employing the differential 
//   cross section by F.W.Lipps & H.A.Tolhoek
//   ( Physica 20 (1954) 395 )
//
#ifndef G4PolarizedComptonCrossSection_h
#define G4PolarizedComptonCrossSection_h 1

#include "G4StokesVector.hh"
#include "G4VPolarizedCrossSection.hh"



class G4PolarizedComptonCrossSection : public G4VPolarizedCrossSection
{
public:
  G4PolarizedComptonCrossSection();
  virtual ~G4PolarizedComptonCrossSection();
public:
  // prepares the ingredients for the calculation of a polarization 
  // dependent differential cross section
  // the kinematics is fixed (X - incoming photon energy in units of electron mass, 
  //  eps - outgoing photon energy in unit of incoming photon energy, 
  // and polarization of the incoming particles fixed (p0, p1)
  // a flag specifies the extent to which polarization is taken 
  // into account
  virtual void Initialize(G4double eps, G4double X, G4double phi,
			  const G4StokesVector & p0,
			  const G4StokesVector & p1,
			  G4int flag=0) override;

  // returns the differential cross section for a given polarisation state
  // of the final state particles to be used in the calculation of the
  // polarization transfer
  // the calculation has to be initialised by calling Initialize()
  // prior to the first call of this function (see above)
  G4double XSection(const G4StokesVector & pol2,const G4StokesVector & pol3) override; 
  // total cross section
  G4double TotalXSection(G4double xmin, G4double xmax, G4double y,
			 const G4StokesVector & pol0,const G4StokesVector & pol1) override;

public:
  // return expected mean polarisation
  G4StokesVector GetPol2() override;
  G4StokesVector GetPol3() override;
private:
  void DefineCoefficients(const G4StokesVector & pol0,
			  const G4StokesVector & pol1);
  // states if an incoming or outgoing particle is polarized
  G4bool gammaPol2, electronPol3;

  // these variables store the information necessary to evaluate the  
  // differential cross section for arbitrary final state 
  // polarizations (used in XSection):
  // - polarization independent part
  G4double phi0;
  // - part depending on the polarization of the final photon
  G4ThreeVector phi2; 
  // - part depending on the polarization of the final electron
  G4ThreeVector phi3;
  // - product of polarizations of initial particles
  G4double polxx, polyy, polzz, polxz, polzx, polyz, polzy, polxy, polyx;

  G4double diffXSFactor, totalXSFactor, re2;
  G4double polXS, unpXS;
};


#endif
