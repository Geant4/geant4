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
// $Id: G4BlineTracerMessenger.hh,v 1.2 2006-06-29 17:15:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
//
// G4BlineTracerMessenger
//
// Class description:
//
// Defines interactive commands that allows to pilot the tracer
// for displaying field lines.

// --------------------------------------------------------------------
// Author: Laurent Desorgher (desorgher@phim.unibe.ch)
//         Created - 2003-10-06
// --------------------------------------------------------------------
#ifndef G4BlineTracerMessenger_h
#define G4BlineTracerMessenger_h 1

#include "G4Types.hh"
#include "G4UImessenger.hh"

class G4BlineTracer;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3Vector;

class G4BlineTracerMessenger : public G4UImessenger
{
  public:  // with description

    G4BlineTracerMessenger(G4BlineTracer* aBlineTool);
    ~G4BlineTracerMessenger();

    void SetNewValue(G4UIcommand * command,G4String newValues);

  private:

    G4BlineTracer* theBlineTool;
    G4UIdirectory* BlineToolDir;

    //  commands

    G4UIcmdWithAnInteger* BlineCmd;
    G4UIcmdWithADoubleAndUnit* SetMaxTrackingStepCmd;
    G4UIcmdWith3Vector* SetDrawColourCmd;
    G4UIcmdWithABool*   SetDrawBlineCmd;    
    G4UIcmdWithABool*  SetDrawPointsCmd;
    G4UIcmdWithADouble* SetPointSizeCmd;
    G4UIcmdWithoutParameter* DrawCmd;
    G4UIcmdWithoutParameter* ResetCmd;
};

#endif
