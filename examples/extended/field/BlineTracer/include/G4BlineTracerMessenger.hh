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
// $Id: G4BlineTracerMessenger.hh,v 1.1 2003/11/25 09:29:46 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
