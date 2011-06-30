//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#ifndef PRIMARYGENERATORMESSENGER_HH
#define PRIMARYGENERATORMESSENGER_HH 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGenerator;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;


class PrimaryGeneratorMessenger : public G4UImessenger {

 public:
   PrimaryGeneratorMessenger(PrimaryGenerator*);
   ~PrimaryGeneratorMessenger();

   void SetNewValue(G4UIcommand*, G4String);

 private:
   PrimaryGenerator* primaryGenerator;

   G4UIdirectory* sourceDirectory;
   G4UIcmdWithADoubleAndUnit* primEnergyCmd;
   G4UIcmdWithADoubleAndUnit* sigmaEnergyCmd;
   G4UIcmdWithADoubleAndUnit* sigmaSpatialCmd;
   G4UIcmdWithADoubleAndUnit* incidAngleCmd;
};

#endif // PRIMARYGENERATORMESSENGER_HH
