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
// $Id: Sc01EventActionMessenger.hh,v 1.2 2006-06-29 18:53:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef Sc01EventActionMessenger_h
#define Sc01EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Sc01EventAction;

class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;



class Sc01EventActionMessenger: public G4UImessenger
{
  public:

    Sc01EventActionMessenger( Sc01EventAction* evAct);
   ~Sc01EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:

    Sc01EventAction*     eventAction;   
    G4UIcmdWithAnInteger* PrintCmd;
};

#endif
