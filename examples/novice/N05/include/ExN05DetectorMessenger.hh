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
// $Id: ExN05DetectorMessenger.hh,v 1.4 2002-01-09 17:24:18 ranjard Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef ExN05DetectorMessenger_h
#define ExN05DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class ExN05DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class ExN05DetectorMessenger: public G4UImessenger
{
  public:
    ExN05DetectorMessenger(ExN05DetectorConstruction * myDet);
    ~ExN05DetectorMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    ExN05DetectorConstruction* myDetector;
    
    G4UIdirectory*             mydetDir;
    G4UIcmdWithAString*        SwitchCmd;
    G4UIcmdWithADoubleAndUnit* TmaxCmd;
    G4UIcmdWithADoubleAndUnit* EminCmd;
    G4UIcmdWithADoubleAndUnit* RminCmd;
};

#endif

