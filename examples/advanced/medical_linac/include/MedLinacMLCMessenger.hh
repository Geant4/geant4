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
// $Id: MedLinacMLCMessenger.hh,v 1.2 2004/11/25 16:28:39 mpiergen Exp $
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
    G4UIcmdWithAString* leafNameCmd;
    G4UIcmdWithADoubleAndUnit* positionCmd;
 };

//****************************************************************************
#endif

