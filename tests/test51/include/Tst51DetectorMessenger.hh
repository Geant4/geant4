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
// $Id: Tst51DetectorMessenger.hh,v 1.1 2005-07-05 11:05:54 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// author: Susanna Guatelli (guatelli@ge.infn.it)
// 
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------

#ifndef Tst51DetectorMessenger_h
#define Tst51DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class Tst51DetectorConstruction;

class Tst51DetectorMessenger: public G4UImessenger
{
public:
  Tst51DetectorMessenger(Tst51DetectorConstruction* );
  ~Tst51DetectorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  Tst51DetectorConstruction* detector;
    
  G4UIdirectory*             tst51Dir;
  G4UIcmdWithAString*        targetMaterialCmd;
  G4UIcmdWithADoubleAndUnit* targetThicknessCmd; 
};
#endif

