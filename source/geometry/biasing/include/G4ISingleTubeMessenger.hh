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
// $Id: G4ISingleTubeMessenger.hh,v 1.2 2002-07-12 10:40:43 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ISingleTubeMessenger
//
// Class description:
//
// a messenger createing the commands to messege a G4ITubeFactory: 
// /imp/cell/<cellname>/zmin <zmin>
// /imp/cell/<cellname>/zmax <zmax>
// /imp/cell/<cellname>/ibase <b>
// /imp/cell/<cellname>/iexpo <e>
// with this commands a cell inside a tube shaped geometry
// can be constructed. The cell will range form zmin to zmax
// and have the importance i = pow(b,e)

//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4ISingleTubeMessenger_hh
#define G4ISingleTubeMessenger_hh G4ISingleTubeMessenger_hh

#include "G4UImessenger.hh"


class G4ITubeFactory;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;

class G4ISingleTubeMessenger : public G4UImessenger {
public:
  G4ISingleTubeMessenger(const G4String &,
			 G4ITubeFactory *);
    // constructor with the G4ITubeFactory to be messeged
   
  void SetNewValue(G4UIcommand * command, G4String newValue);
 
private:
  
  G4String fCellName;
  G4ITubeFactory *fITubeFactory;
  G4bool fICellCreated;

  G4UIcmdWithADoubleAndUnit *fZminCmd;
  G4double fZmin;
  G4bool fZminIsSet;
  G4UIcmdWithADoubleAndUnit *fZmaxCmd;
  G4double fZmax;
  G4bool fZmaxIsSet;
  G4UIcmdWithADouble *fIbaseCmd;
  G4double fIbase;
  G4bool fIbasesIsSet;
  G4UIcmdWithADouble *fIexpoCmd;
  G4double fIexpo;
  G4bool fIexpoIsSet;

};

#endif
