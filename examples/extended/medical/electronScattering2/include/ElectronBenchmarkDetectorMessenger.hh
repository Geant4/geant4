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

#ifndef ELECTRONBENCHMARKDETECTORMESSENGER_HH
#define ELECTRONBENCHMARKDETECTORMESSENGER_HH 1

#include "globals.hh"
#include "G4UImessenger.hh"

class ElectronBenchmarkDetector;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class ElectronBenchmarkDetectorMessenger: public G4UImessenger {

public:
  
  ElectronBenchmarkDetectorMessenger(ElectronBenchmarkDetector* det);
  
  ~ElectronBenchmarkDetectorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  
  ElectronBenchmarkDetector* detector;   
  G4UIdirectory* listDir;
  G4UIcmdWithAString* primFoilMatCmd;
  G4UIcmdWithADoubleAndUnit* primFoilThickCmd;
};

#endif








