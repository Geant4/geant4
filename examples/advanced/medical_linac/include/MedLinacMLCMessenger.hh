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
// $Id: MedLinacMLCMessenger.hh,v 1.1 2004-11-24 16:53:29 mpiergen Exp $
//
//
// Code developed by: M. Piergentili

#ifndef MedLinacMLCMessenger_h
#define MedLinacMLCMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class MedLinacMLCDecorator;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//****************************************************************************

class MedLinacMLCMessenger: public G4UImessenger
{
  public:
    MedLinacMLCMessenger(MedLinacMLCDecorator* );
   ~MedLinacMLCMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    MedLinacMLCDecorator* MedLinacMLCDeco;
    
    G4UIdirectory*             MedLinacMLCDir;
    G4UIdirectory*             aDir;
    G4UIcmdWithADoubleAndUnit* a1yPosCmd;
    G4UIcmdWithADoubleAndUnit* a2yPosCmd;
    G4UIcmdWithADoubleAndUnit* a3yPosCmd;
    G4UIcmdWithADoubleAndUnit* a4yPosCmd;
    G4UIcmdWithADoubleAndUnit* a5yPosCmd;
    G4UIcmdWithADoubleAndUnit* a6yPosCmd;
    G4UIcmdWithADoubleAndUnit* a7yPosCmd;
    G4UIcmdWithADoubleAndUnit* a8yPosCmd;
    G4UIcmdWithADoubleAndUnit* a9yPosCmd;
    G4UIcmdWithADoubleAndUnit* a10yPosCmd;
    G4UIcmdWithADoubleAndUnit* a11yPosCmd;
    G4UIcmdWithADoubleAndUnit* a12yPosCmd;
    G4UIcmdWithADoubleAndUnit* a13yPosCmd;
    G4UIcmdWithADoubleAndUnit* a14yPosCmd;
    G4UIcmdWithADoubleAndUnit* a15yPosCmd;
    G4UIcmdWithADoubleAndUnit* a16yPosCmd;
    G4UIcmdWithADoubleAndUnit* a17yPosCmd;
    G4UIcmdWithADoubleAndUnit* a18yPosCmd;
    G4UIcmdWithADoubleAndUnit* a19yPosCmd;
    G4UIcmdWithADoubleAndUnit* a20yPosCmd;
    G4UIcmdWithADoubleAndUnit* a21yPosCmd;
    G4UIcmdWithADoubleAndUnit* a22yPosCmd;
    G4UIcmdWithADoubleAndUnit* a23yPosCmd;
    G4UIcmdWithADoubleAndUnit* a24yPosCmd;
    G4UIcmdWithADoubleAndUnit* a25yPosCmd;
    G4UIcmdWithADoubleAndUnit* a26yPosCmd;
    G4UIcmdWithADoubleAndUnit* a27yPosCmd;
    G4UIcmdWithADoubleAndUnit* a28yPosCmd;
    G4UIcmdWithADoubleAndUnit* a29yPosCmd;
    G4UIcmdWithADoubleAndUnit* a30yPosCmd;
    G4UIcmdWithADoubleAndUnit* a31yPosCmd;
    G4UIcmdWithADoubleAndUnit* a32yPosCmd;
    G4UIcmdWithADoubleAndUnit* a33yPosCmd;
    G4UIcmdWithADoubleAndUnit* a34yPosCmd;
    G4UIcmdWithADoubleAndUnit* a35yPosCmd;
    G4UIcmdWithADoubleAndUnit* a36yPosCmd;
    G4UIcmdWithADoubleAndUnit* a37yPosCmd;
    G4UIcmdWithADoubleAndUnit* a38yPosCmd;
    G4UIcmdWithADoubleAndUnit* a39yPosCmd;
    G4UIcmdWithADoubleAndUnit* a40yPosCmd;


    G4UIcmdWithADoubleAndUnit* b1yPosCmd;
    G4UIcmdWithADoubleAndUnit* b2yPosCmd;
    G4UIcmdWithADoubleAndUnit* b3yPosCmd;
    G4UIcmdWithADoubleAndUnit* b4yPosCmd;
    G4UIcmdWithADoubleAndUnit* b5yPosCmd;
    G4UIcmdWithADoubleAndUnit* b6yPosCmd;
    G4UIcmdWithADoubleAndUnit* b7yPosCmd;
    G4UIcmdWithADoubleAndUnit* b8yPosCmd;
    G4UIcmdWithADoubleAndUnit* b9yPosCmd;
    G4UIcmdWithADoubleAndUnit* b10yPosCmd;
    G4UIcmdWithADoubleAndUnit* b11yPosCmd;
    G4UIcmdWithADoubleAndUnit* b12yPosCmd;
    G4UIcmdWithADoubleAndUnit* b13yPosCmd;
    G4UIcmdWithADoubleAndUnit* b14yPosCmd;
    G4UIcmdWithADoubleAndUnit* b15yPosCmd;
    G4UIcmdWithADoubleAndUnit* b16yPosCmd;
    G4UIcmdWithADoubleAndUnit* b17yPosCmd;
    G4UIcmdWithADoubleAndUnit* b18yPosCmd;
    G4UIcmdWithADoubleAndUnit* b19yPosCmd;
    G4UIcmdWithADoubleAndUnit* b20yPosCmd;
    G4UIcmdWithADoubleAndUnit* b21yPosCmd;
    G4UIcmdWithADoubleAndUnit* b22yPosCmd;
    G4UIcmdWithADoubleAndUnit* b23yPosCmd;
    G4UIcmdWithADoubleAndUnit* b24yPosCmd;
    G4UIcmdWithADoubleAndUnit* b25yPosCmd;
    G4UIcmdWithADoubleAndUnit* b26yPosCmd;
    G4UIcmdWithADoubleAndUnit* b27yPosCmd;
    G4UIcmdWithADoubleAndUnit* b28yPosCmd;
    G4UIcmdWithADoubleAndUnit* b29yPosCmd;
    G4UIcmdWithADoubleAndUnit* b30yPosCmd;
    G4UIcmdWithADoubleAndUnit* b31yPosCmd;
    G4UIcmdWithADoubleAndUnit* b32yPosCmd;
    G4UIcmdWithADoubleAndUnit* b33yPosCmd;
    G4UIcmdWithADoubleAndUnit* b34yPosCmd;
    G4UIcmdWithADoubleAndUnit* b35yPosCmd;
    G4UIcmdWithADoubleAndUnit* b36yPosCmd;
    G4UIcmdWithADoubleAndUnit* b37yPosCmd;
    G4UIcmdWithADoubleAndUnit* b38yPosCmd;
    G4UIcmdWithADoubleAndUnit* b39yPosCmd;
    G4UIcmdWithADoubleAndUnit* b40yPosCmd;
 };

//****************************************************************************
#endif

