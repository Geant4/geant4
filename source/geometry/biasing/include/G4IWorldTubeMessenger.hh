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
// $Id: G4IWorldTubeMessenger.hh,v 1.2 2002-07-12 10:40:43 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4IWorldTubeMessanger
//
// Class description:
//
// a messenger createing the commands for the specification and 
// construction  of a world tube for a G4ImportanceGeometryConstructor: 
// /imp/worldvolume/tube/radius <r>
// /imp/worldcolume/tube/halfwidth <half_z>
// 
// Messeges a G4ImportanceGeometryConstructor with a G4Tube of radius <r>
// and half width along z of half_z.
// Contruct a factory to create cells inside the logical volume
// which is constructed by G4ImportanceGeometryConstructor according
// to the G4Tube given to it.
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4IWorldTubeMessanger_hh
#define G4IWorldTubeMessanger_hh G4IWorldTubeMessanger_hh

#include "G4UImessenger.hh"

class G4UIcmdWithADoubleAndUnit;
class G4VSolid;
class G4ITubeFactory;
class G4ImportanceGeometryConstructor;
class G4LogicalVolume;
class G4VIStore;

class G4IWorldTubeMessenger : public  G4UImessenger
{
public:
  G4IWorldTubeMessenger(G4ImportanceGeometryConstructor *);
    // constructor with a G4ImportanceGeometryConstructor
    // to be messeged

  void SetNewValue(G4UIcommand * command, G4String newValue);

  void ConstructWorldSolid();
    // construct the specified G4Tube

  void ConstructICellFactory(G4LogicalVolume *,
			     G4VIStore *);
    // construct the cell factory for tubes.

private:

  void Error(const G4String& m){
    G4Exception("Error: G4IWorldTubeMessenger: " + m);
  }
  
  G4ImportanceGeometryConstructor *fIGconst;

  G4VSolid *fWorldSolid;
  G4UIcmdWithADoubleAndUnit *fRadiusCmd;
  G4double fRadius;
  G4bool fRadiusIsSet;
  G4UIcmdWithADoubleAndUnit *fHalfHightCmd;
  G4bool fHalfhightIsSet;
  G4double fHalfHight;
  G4ITubeFactory *fITubeFactory;

};
#endif
