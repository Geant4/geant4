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
// $Id: G4ITubeFactory.hh,v 1.3 2002-08-28 15:16:21 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ITubeFactory
//
// Class description:
//
// This class constructs cells of tube shape inside a mother volume.
// The cells have to be given a name. By that name a importance
// value can be assigned to the cell.
// The class is used to construct cells in a volume
// created by  G4ImportanceGeometryConstructor.
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4ITubeFactory_hh 
#define G4ITubeFactory_hh  G4ITubeFactory_hh 

#include "globals.hh"
#include "g4std/map"

class G4ITubeMessenger;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4IStore;
class G4Material;
typedef G4std::map<G4String, G4Tubs *>  G4MapNameTube;
typedef G4std::map<G4String, G4LogicalVolume *>  G4MapNameLogic;
typedef G4std::map<G4String, G4VPhysicalVolume *>  G4MapNamePhysical;


class G4ITubeFactory {
public:
  G4ITubeFactory(G4LogicalVolume *wl, G4IStore *is,
		 G4double Radius, G4double HalfHight);
    // a factory for cells consisting of slabs with tube shape
    // of radius "Radius" and variable dimension along the z-axis.

  ~G4ITubeFactory();
  void AddCell(const G4String &celname, G4double zmin, G4double zmax);
    // add a cell named <celname> of radius "Radius" ranging from 
    // zmin to zmax.

  void SetImportance(const G4String &celname, G4double b, G4double e);
    // set importance of the cell <celname> as i = pow(b,e)

private:

  void Error(const G4String &m) {
    G4Exception("Error: G4ITubeFactory: " + m);
  }

  G4LogicalVolume *fWorldLogic;
  G4IStore *fIStore;
  G4double fRadius;
  G4double fHalfHight;
  G4ITubeMessenger *fITubeMessenger;
  G4Material *fGalactic;

  G4MapNameTube fMapNameTube;
  G4MapNameLogic fMapNameLogic;
  G4MapNamePhysical fMapNamePhysical;


}; 

#endif
