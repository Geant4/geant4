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
// $Id: OlapSteppingAction.hh,v 1.3 2006-06-29 17:22:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// OlapSteppingAction
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#ifndef OlapSteppingAction_h
#define OlapSteppingAction_h

#include <vector>
#include <iostream>

#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NavigationHistory.hh"

class OlapEventAction;

struct OlapStepInfo
{
   G4ThreeVector thePoint;
   G4NavigationHistory theHist;
   
};

std::ostream& operator << (std::ostream &, const OlapStepInfo &);
std::ostream& operator << (std::ostream &, std::vector<OlapStepInfo*> &);

class OlapSteppingAction : public G4UserSteppingAction
{
public:
   OlapSteppingAction(OlapEventAction*);
   ~OlapSteppingAction();
   
   void UserSteppingAction(const G4Step *);
   
   OlapEventAction * theEventAction;
   
};
#endif
