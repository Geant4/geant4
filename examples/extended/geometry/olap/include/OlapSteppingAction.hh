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
// $Id: OlapSteppingAction.hh,v 1.2 2003/06/16 16:49:22 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
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
