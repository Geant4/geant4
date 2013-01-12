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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 17191/03/NL/LvH (Aurora Programme).                     *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#ifndef G4TripathiLightCrossSection_h
#define G4TripathiLightCrossSection_h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4TripathiLightCrossSection.hh
//
// Version:		B.1
// Date:		15/04/04
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		17191/03/NL/LvH
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 6 October 2003, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// 15 March 2004, P R Truscott, QinetiQ Ltd, UK
// Beta release
//
// 26 Dec 2006, D. Wright - added isotope dependence
//
// 19 Aug 2011 V.Ivanchenko move to new design and make x-section per element
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class Description
// Implementation of formulas of Tripathi, Cucinotta and Wilson, NASA Technical
// Paper TP-1999-209726 to calculate cross-sections for nuclear-nuclear
// inelastic scattering for light nuclear systems.
//
// Class Description - End
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"

#include "G4VCrossSectionDataSet.hh"
#include "G4WilsonRadius.hh"

///////////////////////////////////////////////////////////////////////////////
//
class G4TripathiLightCrossSection : public G4VCrossSectionDataSet
{
public:

  G4TripathiLightCrossSection();
  ~G4TripathiLightCrossSection();

  virtual G4bool IsElementApplicable(const G4DynamicParticle* theProjectile,
				     G4int Z, const G4Material*);

  virtual 
  G4double GetElementCrossSection(const G4DynamicParticle* theProjectile,
				  G4int Z, const G4Material* mat = 0);

  inline void SetLowEnergyCheck(G4bool);

private:

  G4WilsonRadius *theWilsonRadius;
  G4double       r_0;
  G4bool         lowEnergyCheck;
};

inline void 
G4TripathiLightCrossSection::SetLowEnergyCheck (G4bool aLowEnergyCheck)
{
  lowEnergyCheck = aLowEnergyCheck;
}

#endif
