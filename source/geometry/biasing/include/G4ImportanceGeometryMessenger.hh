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
// $Id: G4ImportanceGeometryMessenger.hh,v 1.2 2002-07-12 10:40:43 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ImportanceGeometryMessenger
//
// Class description:
//
// a messenger createing the command: 
// /imp/worldvolume/solid <solidtypename>
// At the moment <solidtypename> may only be tube.
// Messeges a G4ImportanceGeometryMessenger.
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4ImportanceGeometryMessenger_hh
#define G4ImportanceGeometryMessenger_hh G4ImportanceGeometryMessenger_hh

#include "G4UImessenger.hh"

class G4ImportanceGeometryConstructor;
class G4UIcmdWithAString;

class G4ImportanceGeometryMessenger : public  G4UImessenger
{
public:
  G4ImportanceGeometryMessenger(G4ImportanceGeometryConstructor &igeo);
    // the constructor takes an G4ImportanceGeometryConstructor
    // which it messages

  void SetNewValue(G4UIcommand * command, G4String newValue);
private:
  
  G4ImportanceGeometryConstructor &fImpGeoConst;
  G4UIcmdWithAString *fSolidTypeCmd;
  
};
#endif
