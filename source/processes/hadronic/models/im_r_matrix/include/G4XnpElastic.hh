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
//      File name:     G4XNNTotal
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4XNPELASTIC_HH
#define G4XNPELASTIC_HH

#include "globals.hh"
#include "G4CrossSectionPatch.hh"
#include "G4CrossSectionVector.hh"

class G4KineticTrack;

class G4XnpElastic : public G4CrossSectionPatch
{

public:

  G4XnpElastic();

  virtual ~G4XnpElastic();

  G4bool operator==(const G4XnpElastic &right) const;
  G4bool operator!=(const G4XnpElastic &right) const;

  virtual const G4CrossSectionVector* GetComponents() const { return components; } 
 
  virtual G4String Name() const;


protected:


private:  

  G4XnpElastic(const G4XnpElastic &right);
  const G4XnpElastic& operator=(const G4XnpElastic &right);

  G4CrossSectionVector* components;

};

#endif


















