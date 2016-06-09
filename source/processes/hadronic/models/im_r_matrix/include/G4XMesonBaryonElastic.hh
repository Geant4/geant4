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
//      File name:     G4XMesonBaryonElastic
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4XMESONBARYONELASTIC_HH
#define G4XMESONBARYONELASTIC_HH

#include "globals.hh"
#include "G4VCrossSectionSource.hh"
#include "G4CrossSectionVector.hh"

class G4KineticTrack;

class G4XMesonBaryonElastic : public G4VCrossSectionSource
{

public:

  G4XMesonBaryonElastic();

  virtual ~G4XMesonBaryonElastic();

  G4bool operator==(const G4XMesonBaryonElastic &right) const;
  G4bool operator!=(const G4XMesonBaryonElastic &right) const;

  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;
 
  virtual const G4CrossSectionVector* GetComponents() const { return 0; }

  virtual G4bool IsValid(G4double e) const;

  virtual G4String Name() const;


protected:


private:  

  G4XMesonBaryonElastic(const G4XMesonBaryonElastic &right);
  const G4XMesonBaryonElastic& operator=(const G4XMesonBaryonElastic &right);
  
  G4double lowLimit;
  G4double highLimit;
  
};

#endif


















