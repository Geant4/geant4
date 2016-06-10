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
// $Id: G4ScoringMessenger.hh 78057 2013-12-03 08:31:15Z gcosmo $
//
// (HISTORY)
//  03-Sep-2007  T.Aso Command definitions are introduced.
//  01-Nov-2007  M.Asai Class is splitted into two.
//  18-Feb-2008  T.Aso Command for cylindrical scorer.
//  18-Feb-2008  T.Aso Command for color chart display.


#ifndef G4ScoringMessenger_h
#define G4ScoringMessenger_h 1

#include "G4UImessenger.hh"

#include <vector>
#include "G4String.hh"

class G4ScoringManager;
class G4VScoringMesh;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcommand;

typedef std::vector<G4String> G4TokenVec;

// class description:
//
//  This is a concrete class of G4UImessenger which handles the commands for
// G4ScoringManager. 
//

class G4ScoringMessenger: public G4UImessenger
{

  public:
    G4ScoringMessenger(G4ScoringManager * SManager);

    ~G4ScoringMessenger();

    void SetNewValue(G4UIcommand * command,G4String newValues);

    G4String GetCurrentValue(G4UIcommand * command);


  protected:    

    void FillTokenVec(G4String newValues,G4TokenVec& token);

    void MeshBinCommand(G4VScoringMesh* mesh, G4TokenVec& token);

  private:

    G4ScoringManager*        fSMan;
    G4UIdirectory*           scoreDir;
    G4UIcmdWithoutParameter* listCmd;
    G4UIcmdWithoutParameter* dumpCmd;
    G4UIcmdWithAnInteger*    verboseCmd;
    //
    G4UIdirectory*             meshCreateDir;
    G4UIcmdWithAString*        meshBoxCreateCmd;
   G4UIcmdWithAString*        meshCylinderCreateCmd;
//    G4UIcmdWithAString*        meshSphereCreateCmd;
    //
    // Mesh commands
    G4UIdirectory*             meshDir;
    G4UIcmdWithAString*        meshOpnCmd;
    //
    G4UIcmdWithoutParameter*   meshClsCmd;
//    G4UIcmdWithABool*        meshActCmd;
    //
    //   Size commands
    G4UIcmdWith3VectorAndUnit* mBoxSizeCmd;
    G4UIcommand* mCylinderSizeCmd;
//    G4UIcmdWith3VectorAndUnit* mSphereSizeCmd;
    //
    //   Division command
    G4UIcommand*               mBinCmd;
    //
    //   Placement command
    G4UIdirectory*             mTransDir;
    G4UIcmdWithoutParameter*   mTResetCmd;
    G4UIcmdWith3VectorAndUnit* mTXyzCmd;
    G4UIdirectory*             mRotDir;
    G4UIcmdWithoutParameter*   mRResetCmd;
    G4UIcmdWithADoubleAndUnit* mRotXCmd;
    G4UIcmdWithADoubleAndUnit* mRotYCmd;
    G4UIcmdWithADoubleAndUnit* mRotZCmd;
    //
    // Draw Command
    G4UIcommand *             drawCmd;
    G4UIcommand *             drawColumnCmd;
    G4UIdirectory *           colorMapDir;
    G4UIcmdWithoutParameter * listColorMapCmd;
    G4UIcmdWithAString *      floatMinMaxCmd;
    G4UIcommand *             colorMapMinMaxCmd;
    //G4UIcommand *             chartCmd;
    //
    // Dump scoring result to file
    G4UIcommand * dumpQtyToFileCmd;
    G4UIcommand * dumpAllQtsToFileCmd;
    //

};




#endif

