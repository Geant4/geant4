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
//
// $Id: G4VisCommandsGeometry.hh 66870 2013-01-14 23:38:59Z adotti $

// /vis/geometry commands - John Allison  31st January 2006

#ifndef G4VISCOMMANDSGEOMETRY_HH
#define G4VISCOMMANDSGEOMETRY_HH

#include "G4VVisCommand.hh"

class G4UIcmdWithAString;
class G4LogicalVolume;
class G4VisAttributes;

#include <map>

class G4VVisCommandGeometry: public G4VVisCommand {
public:
  virtual ~G4VVisCommandGeometry();
protected:
  static std::map<G4LogicalVolume*, const G4VisAttributes*>
  fVisAttsMap;
  typedef std::map<G4LogicalVolume*, const G4VisAttributes*>::const_iterator
  VisAttsMapIterator;
};

class G4VisCommandGeometryList: public G4VVisCommandGeometry {
public:
  G4VisCommandGeometryList();
  virtual ~G4VisCommandGeometryList();
  G4String GetCurrentValue(G4UIcommand* command);
  void SetNewValue(G4UIcommand* command, G4String newValue);
private:
  G4VisCommandGeometryList(const G4VisCommandGeometryList&);
  G4VisCommandGeometryList& operator=(const G4VisCommandGeometryList&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandGeometryRestore: public G4VVisCommandGeometry {
public:
  G4VisCommandGeometryRestore();
  virtual ~G4VisCommandGeometryRestore();
  G4String GetCurrentValue(G4UIcommand* command);
  void SetNewValue(G4UIcommand* command, G4String newValue);
private:
  G4VisCommandGeometryRestore(const G4VisCommandGeometryRestore&);
  G4VisCommandGeometryRestore& operator=(const G4VisCommandGeometryRestore&);
  G4UIcmdWithAString* fpCommand;
};

#endif
