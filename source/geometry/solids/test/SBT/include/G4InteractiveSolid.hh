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
class G4UIcmdWithAnInteger;

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

   // Constituent Solids for Boolean
  G4VSolid	*fSolidA;
  G4VSolid	*fSolidB;
	
  G4UIdirectory	*volumeDirectory;
	
  G4UIcmdParg		*boxArgs[3];
  G4UIcmdWithPargs	*boxCmd;
  void MakeMeABox( G4String values );
	
  G4UIcmdParg		*consArgs[7];
  G4UIcmdWithPargs	*consCmd;
  void MakeMeACons( G4String values );
	
  G4UIcmdParg		*orbArgs[1];
  G4UIcmdWithPargs	*orbCmd;
  void MakeMeAnOrb( G4String values );
	
  G4UIcmdParg		*paraArgs[6];
  G4UIcmdWithPargs	*paraCmd;
  void MakeMeAPara( G4String values );
	
  G4UIcmdParg		*sphereArgs[6];
  G4UIcmdWithPargs	*sphereCmd;
  void MakeMeASphere( G4String values );
	
  G4UIcmdParg		*torusArgs[5];
  G4UIcmdWithPargs	*torusCmd;
  void MakeMeATorus( G4String values );
	
  G4UIcmdParg		*trapArgs[11];
  G4UIcmdWithPargs	*trapCmd;
  void MakeMeATrap( G4String values );
	
  G4UIcmdParg		*trdArgs[5];
  G4UIcmdWithPargs	*trdCmd;
  void MakeMeATrd( G4String values );
 
  G4UIcmdParg		*gentrapArgs[3];
  G4UIcmdWithPargs	*gentrapCmd;
  void MakeMeAGenericTrap( G4String values );

  G4UIcmdParg		*parabolArgs[3];
  G4UIcmdWithPargs	*parabolCmd;
  void MakeMeAParaboloid( G4String values );
	
  G4UIcmdParg		*tubsArgs[5];
  G4UIcmdWithPargs	*tubsCmd;
  void MakeMeATubs( G4String values );

  G4UIcmdParg		*cuttubsArgs[7];
  G4UIcmdWithPargs	*cuttubsCmd;
  void MakeMeACutTubs( G4String values );
	
  G4UIcmdParg		*ellipsoidArgs[5];
  G4UIcmdWithPargs	*ellipsoidCmd;
  void MakeMeAnEllipsoid( G4String values );
	
  G4UIcmdParg		*elConeArgs[4];
  G4UIcmdWithPargs	*elConeCmd;
  void MakeMeAnEllipticalCone( G4String values );
	
  G4UIcmdParg		*elTubeArgs[6];
  G4UIcmdWithPargs	*elTubeCmd;
  void MakeMeAnEllipticalTube( G4String values );

  G4UIcmdParg		*extrudedArgs[8];
  G4UIcmdWithPargs	*extrudedCmd;
  void MakeMeAnExtrudedSolid( G4String values );

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
	
  G4UIcmdParg		*tesselArgs[9];
  G4UIcmdWithPargs	*tesselCmd;
  void MakeMeATessellatedSolid( G4String values );
	
  G4UIcmdParg		*tessel2Args[8];
  G4UIcmdWithPargs	*tessel2Cmd;
  void MakeMeATessellatedSolid2( G4String values );

  G4UIcmdParg		*tetArgs[4];
  G4UIcmdWithPargs	*tetCmd;
  void MakeMeATet( G4String values );
	
  G4UIcmdParg		*twistedBoxArgs[4];
  G4UIcmdWithPargs	*twistedBoxCmd;
  void MakeMeATwistedBox( G4String values );
	
  G4UIcmdParg		*twistedTrapArgs[5];
  G4UIcmdWithPargs	*twistedTrapCmd;
  void MakeMeATwistedTrap( G4String values );
	
  G4UIcmdParg		*twistedTrap2Args[11];
  G4UIcmdWithPargs	*twistedTrap2Cmd;
  void MakeMeATwistedTrap2( G4String values );
	
  G4UIcmdParg		*twistedTrdArgs[6];
  G4UIcmdWithPargs	*twistedTrdCmd;
  void MakeMeATwistedTrd( G4String values );
	
  G4UIcmdParg		*twistedTubsArgs[7];
  G4UIcmdWithPargs	*twistedTubsCmd;
  void MakeMeATwistedTubs( G4String values );
	
  G4UIcmdWithPargs	*dircTestCmd;
  void MakeMeDircTest();

  typedef enum BooleanOp {
    INTERSECTION,
    SUBTRACTION,
    UNION
  } BooleanOp;

  // G4UIcmdWithPargs	*TestBooleanSolidCmd;  
  // void MakeMeTestBooleanSolid(G4String values);

/* 
  G4UIcmdWithAnInteger  *StoreConstituentSolidCmd;
  void StoreConstituentSolid(G4String values);
*/

  // G4UIcmdParg		*SimpleBooleanArgs[1];
  // G4UIcmdPargInteger    *SimpleBooleanType;
  G4UIcmdWithAnInteger  *SimpleBooleanSolidCmd;
  void MakeMeASimpleBooleanSolid(G4String values);

  // G4UIcmdWithPargs      *GeneralBooleanSolidCmd;
  // void MakeMeGeneralBooleanSolid(G4String values);

  /* Here add new commands and functions to create solids */
};


#endif
