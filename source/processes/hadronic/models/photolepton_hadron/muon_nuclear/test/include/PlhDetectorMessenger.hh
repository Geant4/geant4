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
// $Id: PlhDetectorMessenger.hh,v 1.1 2003-11-12 17:21:55 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef PlhDetectorMessenger_h
#define PlhDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PlhDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class PlhDetectorMessenger: public G4UImessenger
{
  public:
    PlhDetectorMessenger(PlhDetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    PlhDetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selMatCmd;

};

#endif

