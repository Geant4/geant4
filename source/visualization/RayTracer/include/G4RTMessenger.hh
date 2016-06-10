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
// $Id: G4RTMessenger.hh 74050 2013-09-20 09:38:19Z gcosmo $
//
//

// class description:
//
//  This is a concrete class of G4UImessenger. This class defines various
// UI commands which are unique to G4RayTracer.
//


#ifndef G4RTMessenger_HH
#define G4RTMessenger_HH 1

#include "G4UImessenger.hh"
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4TheRayTracer;

class G4RTMessenger : public G4UImessenger
{
  public:
    static G4RTMessenger* GetInstance(G4TheRayTracer* p1);  // Singleton constructor.
    virtual ~G4RTMessenger();
    
    virtual G4String GetCurrentValue(G4UIcommand * command);
    virtual void SetNewValue(G4UIcommand * command,G4String newValue);

  private:
    G4RTMessenger(G4TheRayTracer* p1);
    static G4RTMessenger* fpInstance;
    G4TheRayTracer* theDefaultTracer;  // The first tracer to
				       // instantiate this messenger.
    G4TheRayTracer* theTracer;         // The current tracer.

    G4UIdirectory* rayDirectory;
    G4UIcmdWithAnInteger* columnCmd;
    G4UIcmdWithAnInteger* rowCmd;
    G4UIcmdWith3VectorAndUnit* targetCmd;
    G4UIcmdWith3VectorAndUnit* eyePosCmd;
    G4UIcmdWith3Vector* lightCmd;
    G4UIcmdWithADoubleAndUnit* spanXCmd;
    G4UIcmdWithADoubleAndUnit* headCmd;
    G4UIcmdWithADoubleAndUnit* attCmd;
    G4UIcmdWithABool* distCmd;
    G4UIcmdWithABool* transCmd;
    G4UIcmdWithAString* fileCmd;
    G4UIcmdWith3Vector* bkgColCmd;
};

#endif



