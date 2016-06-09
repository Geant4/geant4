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
//    *****************************************
//    *                                       *
//    *      BrachyDetectrorMessenger.hh      *
//    *                                       *
//    *****************************************
//
// $Id: BrachyDetectorMessenger.hh,v 1.5 2003/05/22 17:20:41 guatelli Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
#ifndef BrachyDetectorMessenger_h
#define BrachyDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class BrachyDetectorConstruction;
class BrachyFactoryIr;
class BrachyRunAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class BrachyDetectorMessenger: public G4UImessenger
{
public:
  BrachyDetectorMessenger(BrachyDetectorConstruction* );
  ~BrachyDetectorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
  
private:

  BrachyDetectorConstruction*  detector;//pointer to detector
  G4UIdirectory*               detectorDir; 
  G4UIcmdWithAString*          phantomMaterialCmd; // change phantom material
  G4UIcmdWithAString*          sourceCmd; //change brachytherapis source 
};
#endif

