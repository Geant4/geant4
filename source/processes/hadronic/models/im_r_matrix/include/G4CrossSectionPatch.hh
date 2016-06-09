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
//
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4CrossSectionPatch
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4CROSSSECTIONPATCH_HH
#define G4CROSSSECTIONPATCH_HH

#include "globals.hh"
#include "G4CrossSectionVector.hh"
#include "G4VCrossSectionSource.hh"

class G4KineticTrack;

class G4CrossSectionPatch : public G4VCrossSectionSource
{

public:

  G4CrossSectionPatch();

  virtual ~G4CrossSectionPatch();

  G4bool operator==(const G4CrossSectionPatch& right) const;
  G4bool operator!=(const G4CrossSectionPatch& right) const;

  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;
 
  virtual const G4CrossSectionVector* GetComponents() const = 0;

  virtual G4bool IsValid(G4double e) const;


protected:

  // Transition between two cross section formulae in an energy range
  G4double Transition(const G4KineticTrack& trk1, 
		      const G4KineticTrack& trk2,
		      const G4VCrossSectionSource* comp1, 
		      const G4VCrossSectionSource* comp2) const;

  G4double Transition(G4double ecm, 
		      G4double sigma1, G4double sigma2, 
		      G4double e1, G4double e2) const;

private:  

  G4CrossSectionPatch(const G4CrossSectionPatch& right);

  G4CrossSectionPatch& operator=(const G4CrossSectionPatch& right);
};
  
#endif
