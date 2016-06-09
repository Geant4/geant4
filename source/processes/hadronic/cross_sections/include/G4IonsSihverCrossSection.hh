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
#ifndef G4IonsSihverCrossSection_h
#define G4IonsSihverCrossSection_h
//
// Class Description
// Implementation of formulas 
// Sihver et al. Phys. Rev. C 47 1225 (1993); 
// Total Reaction Cross Section for Nucleus-nucles reactions.
//    Energy independent   
//    Valid for 100>MeV/nucleon 
// Class Description - End
//
// 23-Dec-2006 Isotope dependence added by D. Wright
// 19-Aug-2011 V.Ivanchenko move to new design and make x-section per element
//

#include "globals.hh"
#include "G4Proton.hh"

#include "G4VCrossSectionDataSet.hh"

class G4IonsSihverCrossSection : public G4VCrossSectionDataSet
{
public:

  G4IonsSihverCrossSection();

  virtual ~G4IonsSihverCrossSection();
   
  virtual
  G4bool IsElementApplicable(const G4DynamicParticle* aDP, G4int Z,
			     const G4Material*);

  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, 
				  G4int Z, const G4Material*);

  virtual void CrossSectionDescription(std::ostream&) const;

private:
  const G4double square_r0;

};

#endif
