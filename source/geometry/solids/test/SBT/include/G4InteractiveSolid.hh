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
// G4InteractiveSolid.hh
//
// A messenger that allows one to construct a solid interactively
// (i.e. via command line)
//
// The solid thus created can be recovered using the G4SolidQuery
// method GetSolid.
//
// Notes:
//    * The G4UIcommand design is somewhat inflexible. It would have
//      been much better to specify each command argument as a
//      class of it's own (with input methods, help methods, etc.)
//      then to create the monolithic monstrosities like G4UICmdWith****.
//      Alas, I'm tempted to fix this, but for the sake of expediency
//      I will cheat and use my own string interpretator. The 
//      side effect is that interactive help is much compromised.
//

#ifndef G4InteractiveSolid_hh
#define G4InteractiveSolid_hh

#include "G4UImessenger.hh"
#include "G4SolidQuery.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithPargs;
class G4UIcmdParg;

class G4InteractiveSolid : public G4UImessenger, public G4SolidQuery {
public:
  G4InteractiveSolid( const G4String &commandPrefix );
  virtual ~G4InteractiveSolid();
	
  inline G4VSolid *GetSolid() const { return solid; }
	
  void SetNewValue( G4UIcommand *command, G4String newValues );
  G4String GetCurrentValue( G4UIcommand *command );
	
protected:
  void DeleteArgArray( G4UIcmdParg **array, const G4int nItem );
  G4String ConvertArgsToString( G4UIcmdParg **array, const G4int nItem );
	
  G4VSolid	*solid;
	
  G4UIdirectory	*volumeDirectory;
	
  G4UIcmdParg		*boxArgs[3];
  G4UIcmdWithPargs	*boxCmd;
  void MakeMeABox( G4String values );
	
  G4UIcmdParg		*paraArgs[6];
  G4UIcmdWithPargs	*paraCmd;
  void MakeMeAPara( G4String values );
	
  G4UIcmdParg		*trapArgs[11];
  G4UIcmdWithPargs	*trapCmd;
  void MakeMeATrap( G4String values );
	
  G4UIcmdParg		*trdArgs[5];
  G4UIcmdWithPargs	*trdCmd;
  void MakeMeATrd( G4String values );
	
  G4UIcmdParg		*sphereArgs[6];
  G4UIcmdWithPargs	*sphereCmd;
  void MakeMeASphere( G4String values );
	
  G4UIcmdParg		*torusArgs[5];
  G4UIcmdWithPargs	*torusCmd;
  void MakeMeATorus( G4String values );
	
  G4UIcmdParg		*tubsArgs[5];
  G4UIcmdWithPargs	*tubsCmd;
  void MakeMeATubs( G4String values );
	
  G4UIcmdParg		*consArgs[7];
  G4UIcmdWithPargs	*consCmd;
  void MakeMeACons( G4String values );
	
  G4UIcmdParg		*hypeArgs[5];
  G4UIcmdWithPargs	*hypeCmd;
  void MakeMeAHype( G4String values );
	
  G4UIcmdParg		*polyconeArgs[5];
  G4UIcmdWithPargs	*polyconeCmd;
  void MakeMeAPolycone( G4String values );
	
  G4UIcmdParg		*polycone2Args[6];
  G4UIcmdWithPargs	*polycone2Cmd;
  void MakeMeAPolycone2( G4String values );
	
  G4UIcmdParg		*polyhedraArgs[6];
  G4UIcmdWithPargs	*polyhedraCmd;
  void MakeMeAPolyhedra( G4String values );
	
  G4UIcmdParg		*polyhedra2Args[7];
  G4UIcmdWithPargs	*polyhedra2Cmd;
  void MakeMeAPolyhedra2( G4String values );
	
  G4UIcmdParg		*ellipticalTubeArgs[6];
  G4UIcmdWithPargs	*ellipticalTubeCmd;
  void MakeMeAnEllipticalTube( G4String values );
	
  G4UIcmdWithPargs	*dircTestCmd;
  void MakeMeDircTest();


  
  typedef enum BooleanOp {
    INTERSECTION,
    SUBTRACTION,
    UNION
  } BooleanOp;

  G4UIcmdWithPargs	*BooleanSolid1Cmd;
  void MakeMeBooleanSolid1(G4String values);

  /* Here add new commands and functions to create solids */

};


#endif
