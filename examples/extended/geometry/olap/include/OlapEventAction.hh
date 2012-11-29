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
/// \file geometry/olap/include/OlapEventAction.hh
/// \brief Definition of the OlapEventAction class
//
//
// $Id$
//
// 
// --------------------------------------------------------------
// OlapEventAction
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#ifndef OlapEventAction_h
#define OlapEventAction_h

#include <vector>
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4NavigationHistory.hh"
#include "OlapSteppingAction.hh"

class OlapRunAction;
class OlapDetConstr;

// class to collect detected overlaps
class OlapInfo
{
public:
   OlapInfo(G4NavigationHistory & h1,
            G4NavigationHistory & h2,
            G4ThreeVector & p1,
            G4ThreeVector & p2,
            G4int a=0,
            G4LogicalVolume* original=0)
     : hist1(h1), hist2(h2), v1(p1), v2(p2),axis(a),
       info(""), probNot(false), originalMother(original) {}
   
   ~OlapInfo();
   
   G4bool operator==(const OlapInfo &);
   
   G4NavigationHistory hist1, hist2;
   G4ThreeVector v1, v2;
   G4int axis;
     // geantino ray travelling parallel to that axis detected the olap
   std::vector<OlapStepInfo *> stAB;
   std::vector<OlapStepInfo *> stBA;
   G4String info;
   G4bool probNot ;            
   G4LogicalVolume* originalMother;
};

std::ostream & operator<<(std::ostream& flux,  OlapInfo & oi);

// ------------==============-------------==============----------------

class OlapEventAction : public G4UserEventAction
{
public:
   OlapEventAction(OlapRunAction*);
   ~OlapEventAction();
   
   void BeginOfEventAction(const G4Event* anEvent);
   void EndOfEventAction(const G4Event* anEvent);
 
   // returns true, if new olverlap was detected!
   G4bool InsertOI(OlapInfo *);
   
   OlapRunAction * theRunAction; 
   OlapDetConstr * theDet;
   std::vector<OlapStepInfo*> ABSteps;
   std::vector<OlapStepInfo*> BASteps;
   G4bool dontDelete;
};

#endif
