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
// $Id: G4PolarizedMollerCrossSection.hh,v 1.2 2006-11-17 14:14:19 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// GEANT4 Class file
//  
//
// File name:     G4PolarizedMollerCrossSection
//
// Author:        Andreas Schaelicke
//
// Creation date: 12.01.2006
//
// Modifications:
//   16-01-06 included cross section as calculated by P.Starovoitov
//
// Class Description:
//   * calculates the differential cross section
//     incomming electron (along positive z direction) scatters at an electron at rest
//   * phi denotes the angle between the scattering plane (defined by the
//     outgoing electron) and X-axis
//   * all stokes vectors refer to spins in the Global System (X,Y,Z)
//
#ifndef G4PolarizedMollerCrossSection_h
#define G4PolarizedMollerCrossSection_h 1

#include "G4StokesVector.hh"
#include "G4VPolarizedCrossSection.hh"

class G4PolarizedMollerCrossSection : public G4VPolarizedCrossSection
{
public:
  G4PolarizedMollerCrossSection();
  virtual ~G4PolarizedMollerCrossSection();
public:
  void Initialize(G4double x, G4double y, G4double phi,
		  const G4StokesVector & p0,const G4StokesVector & p1,
		  G4int flag=0);

  G4double XSection(const G4StokesVector & pol2,const G4StokesVector & pol3);
  G4double TotalXSection(G4double xmin, G4double xmax, G4double y,
			 const G4StokesVector & pol0,const G4StokesVector & pol1);
  // return expected mean polarisation
  G4StokesVector GetPol2();
  G4StokesVector GetPol3();
private:
  G4double phi0;
  // - part depending on the polarization of the final electron P1
  G4ThreeVector phi2;
  // - part depending on the polarization of the final electron P2
  G4ThreeVector phi3;

};
#endif
