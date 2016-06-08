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
//      File name:     G4CrossSectionComposite
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4CROSSSECTIONCOMPOSITE_HH
#define G4CROSSSECTIONCOMPOSITE_HH

#include "globals.hh"
#include "G4CrossSectionVector.hh"
#include "G4VCrossSectionSource.hh"

class G4KineticTrack;

class G4CrossSectionComposite : public G4VCrossSectionSource
{

public:

  G4CrossSectionComposite();

  virtual ~G4CrossSectionComposite();

  G4bool operator==(const G4CrossSectionComposite& right) const;
  G4bool operator!=(const G4CrossSectionComposite& right) const;

  // Cross section of composite is the sum of components cross sections
  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;
 
  virtual const G4CrossSectionVector* GetComponents() const = 0;

  virtual G4bool IsValid(G4double e) const;


protected:

private:  

  G4CrossSectionComposite(const G4CrossSectionComposite& right);

  G4CrossSectionComposite& operator=(const G4CrossSectionComposite& right);
};
  
#endif
