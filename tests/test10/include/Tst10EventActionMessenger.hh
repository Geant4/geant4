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
// $Id: Tst10EventActionMessenger.hh,v 1.1 2004-01-25 14:04:28 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef Tst10EventActionMessenger_h
#define Tst10EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst10EventAction;

class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;



class Tst10EventActionMessenger: public G4UImessenger
{
  public:

    Tst10EventActionMessenger( Tst10EventAction* evAct);
   ~Tst10EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:

    Tst10EventAction*     eventAction;   
    G4UIcmdWithAnInteger* PrintCmd;
};

#endif
