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
// $Id: G4GeometryMessenger.hh,v 1.2 2001-10-24 18:10:38 gcosmo Exp $
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

#ifndef G4GeometryMessenger_hh
#define G4GeometryMessenger_hh

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ThreeVector.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4TransportationManager;

class G4GeometryMessenger : public G4UImessenger
{
  public:  // with description

    G4GeometryMessenger(G4TransportationManager* tman);
    ~G4GeometryMessenger();
      // Constructor and destructor

    void SetNewValue( G4UIcommand* command, G4String newValues );
    G4String GetCurrentValue( G4UIcommand* command );
  
  protected:

    void ResetNavigator();
    void RunTest();

    G4UIdirectory             *geodir, *navdir, *testdir;
    G4UIcmdWith3VectorAndUnit *posCmd, *dirCmd;
    G4UIcmdWithoutParameter   *runCmd, *resCmd;
    G4UIcmdWithABool          *addCmd;
    G4UIcmdWithADoubleAndUnit *tolCmd;
  
  private:

    G4bool geometryOpened;
    G4ThreeVector x, p;
    G4bool        extra, newtol;
    G4double      tol;

  private:
  
    G4TransportationManager* tmanager;

};

#endif
