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
// $Id: G4VPolarizedCrossSection.hh,v 1.3 2007/07/10 09:36:39 schaelic Exp $
// GEANT4 tag $Name: geant4-09-00-patch-02 $
//
// File name:     G4VPolarizedCrossSection
//
// Author:        Andreas Schaelicke
//
// Creation date: 15.05.2005
//
// Modifications:
//
// Class Description:
//   (pure virtual) interface class
//
//   provides readable but efficient routines to determine 
//   polarization for the final state of a given process
//   empoying the differential cross section
//

#ifndef G4VPolarizedCrossSection_h
#define G4VPolarizedCrossSection_h 1

#include "G4StokesVector.hh"


class G4VPolarizedCrossSection
{
public:
  G4VPolarizedCrossSection();
  virtual ~G4VPolarizedCrossSection();
public:
  virtual void Initialize(G4double, G4double, G4double,
			  const G4StokesVector & p0,const G4StokesVector & p1,
 			  G4int flag=0); 
  virtual G4double XSection(const G4StokesVector & pol2,const G4StokesVector & pol3) = 0; 
  virtual G4double TotalXSection(G4double xmin, G4double xmax, G4double y,
				 const G4StokesVector & pol0,
				 const G4StokesVector & pol1);
  
  // return expected mean polarisation
  virtual G4StokesVector GetPol2();
  virtual G4StokesVector GetPol3();

  // return appropriate distribute polarisation states;
//  void DicePolarization();
//   G4StokesVector DicedPol2();
//   G4StokesVector DicedPol3();
protected:
//   G4double costheta; 
//  G4StokesVector pol0, pol1;
//   G4StokesVector pol2, pol3;
};


#endif
