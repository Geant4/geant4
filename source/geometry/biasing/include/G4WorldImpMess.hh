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
// $Id: G4WorldImpMess.hh,v 1.2 2002-07-12 10:40:43 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4WorldImpMess
//
// Class description:
//
// a messenger createing the commands: 
// /imp/worldvolume/ibase <b>
// /imp/worldvolume/iexpo <e>
// 
// Messeges a G4ImportanceGeometryConstructor with the base <b> and 
// exponent <e> of the importance value for the world volume.
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4WorldImpMess_hh
#define G4WorldImpMess_hh G4WorldImpMess_hh

#include "G4UImessenger.hh"

class G4ImportanceGeometryConstructor;
class G4UIcmdWithADouble;

class G4WorldImpMess : public G4UImessenger {
public:
  G4WorldImpMess(G4ImportanceGeometryConstructor *);
    // constructed with the G4ImportanceGeometryConstructor
    // to be messeged.
 
  void SetNewValue(G4UIcommand * command, G4String newValue);
private:
  G4ImportanceGeometryConstructor *fIGConst;

  G4UIcmdWithADouble *fIbaseCmd;
  G4double fIbase;
  G4bool fIbaseIsSet;
  G4UIcmdWithADouble *fIexpoCmd;
  G4double fIexpo;
  G4bool fIexpoIsSet;

};

#endif
