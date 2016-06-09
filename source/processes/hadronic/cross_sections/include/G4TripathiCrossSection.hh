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
// by JPW, working, but to be cleaned up. @@@@

#ifndef G4TripathiCrossSection_h
#define G4TripathiCrossSection_h

// 19-Aug-2011 V.Ivanchenko move to new design and make x-section per element

// Class Description
// Implementation of formulas in analogy to NASA technical paper 3621 by 
// Tripathi, et al.; Cross-sections for ion ion scattering.

// Class Description - End

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"

class G4TripathiCrossSection : public G4VCrossSectionDataSet
{
public:
   
  G4TripathiCrossSection();
  ~G4TripathiCrossSection();

  virtual
  G4bool IsElementApplicable(const G4DynamicParticle* aPart, 
			     G4int Z, const G4Material*);

  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, 
				  G4int Z, const G4Material*);

};

#endif
