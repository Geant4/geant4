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
// /vis/scene/plotter commands - Guy Barrand September 2021.
//

#ifndef G4VISCOMMANDSPLOTTER_HH
#define G4VISCOMMANDSPLOTTER_HH

#include "G4VVisCommand.hh"

#define G4VIS_COMMAND_PLOTTER_HEADER(a__cmd)\
class G4VisCommandPlotter##a__cmd: public G4VVisCommand {\
public:\
  G4VisCommandPlotter##a__cmd ();\
  virtual ~G4VisCommandPlotter##a__cmd ();\
  G4String GetCurrentValue (G4UIcommand*) {return "";}\
  void SetNewValue (G4UIcommand* command, G4String newValue);\
private:\
  G4VisCommandPlotter##a__cmd (const G4VisCommandPlotter##a__cmd&);\
  G4VisCommandPlotter##a__cmd& operator = (const G4VisCommandPlotter##a__cmd&);\
\
  G4UIcommand* fpCommand;\
};

G4VIS_COMMAND_PLOTTER_HEADER(Create)
G4VIS_COMMAND_PLOTTER_HEADER(SetLayout)
G4VIS_COMMAND_PLOTTER_HEADER(AddStyle)
G4VIS_COMMAND_PLOTTER_HEADER(AddRegionStyle)
G4VIS_COMMAND_PLOTTER_HEADER(AddRegionParameter)
G4VIS_COMMAND_PLOTTER_HEADER(Clear)
G4VIS_COMMAND_PLOTTER_HEADER(ClearRegion)
G4VIS_COMMAND_PLOTTER_HEADER(List)

G4VIS_COMMAND_PLOTTER_HEADER(AddRegionH1)
G4VIS_COMMAND_PLOTTER_HEADER(AddRegionH2)

#undef G4VIS_COMMAND_PLOTTER_HEADER

#endif
