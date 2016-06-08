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
