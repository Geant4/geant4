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
// R&D: Vladimir.Grichine@cern.ch
//      Simone.Gilardoni@cern.ch


#ifndef G4HadFileSpec_HH
#define G4HadFileSpec_HH 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"

class G4HadFileSpec
{
public:

  G4HadFileSpec();
  G4HadFileSpec(G4String&, G4Isotope*, G4String&,G4String&);
  G4HadFileSpec(G4String&, G4Element*, G4String&,G4String&);
  G4HadFileSpec(G4String&, G4Material*, G4String&,G4String&);


  G4String G4HDSFilename();
  G4String G4HDSFilepath();


private:
  
  G4String&  fprimary;
  G4Isotope* fisotope;
  G4Element*  felement;
  G4Material* fmaterial;
  G4String&  fsecondary;
  G4String&  fprocess;
 
  
};
#endif
