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
// $Id: G4GeometryMessenger.hh,v 1.5 2010-11-10 14:06:40 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4GeometryMessenger
//
// Class description:
//
// A messenger defining commands for debugging, verifying
// and controlling the detector geometry and navigation.

// Author: G.Cosmo, CERN.
// --------------------------------------------------------------------
#ifndef G4GeometryMessenger_hh
#define G4GeometryMessenger_hh

#include "G4Types.hh"
#include "G4UImessenger.hh"
#include "G4ThreeVector.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4TransportationManager;
class G4GeomTestStreamLogger;
class G4GeomTestVolume;

class G4GeometryMessenger : public G4UImessenger
{
  public:  // with description

    G4GeometryMessenger(G4TransportationManager* tman);
    ~G4GeometryMessenger();
      // Constructor and destructor

    void SetNewValue( G4UIcommand* command, G4String newValues );
    G4String GetCurrentValue( G4UIcommand* command );
  
  private:

    void Init();
    void CheckGeometry();
    void ResetNavigator();
    void SetVerbosity(G4String newValue);
    void SetCheckMode(G4String newValue);
    void SetPushFlag(G4String newValue);
    void LineTest();
    void RecursiveLineTest();
    void GridTest();
    void RecursiveGridTest();
    void CylinderTest();
    void RecursiveCylinderTest();

    G4UIdirectory             *geodir, *navdir, *testdir;
    G4UIcmdWith3VectorAndUnit *posCmd, *dirCmd;
    G4UIcmdWith3Vector        *grzCmd, *cyzCmd;
    G4UIcmdWithABool          *chkCmd, *pchkCmd, *linCmd,
                              *grdCmd, *cylCmd, *runCmd;
    G4UIcmdWithoutParameter   *recCmd, *resCmd;
    G4UIcmdWithADoubleAndUnit *tolCmd;
    G4UIcmdWithAnInteger      *verbCmd, *rcsCmd, *rcdCmd;
    G4UIcmdWithADouble        *cfzCmd, *cfrCmd;

    G4ThreeVector x, p, grdRes, cylRes;
    G4double      cylfZ, cylfR;
    G4bool        newtol;
    G4double      tol;
    G4int         recLevel, recDepth;

    G4TransportationManager* tmanager;
    G4GeomTestStreamLogger* tlogger;
    G4GeomTestVolume* tvolume;
};

#endif
