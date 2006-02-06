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
// $Id: G4VisCommandsGeometry.hh,v 1.1 2006-02-06 12:13:52 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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
  static std::map<const G4LogicalVolume*, const G4VisAttributes*>
  fVisAttsMap;
  typedef
  std::map<const G4LogicalVolume*, const G4VisAttributes*>::const_iterator
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
