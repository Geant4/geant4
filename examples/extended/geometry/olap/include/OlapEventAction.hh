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
// $Id: OlapEventAction.hh,v 1.2 2003/06/16 16:49:17 gunter Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
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
	    G4LogicalVolume* original=0) :
	    hist1(h1), hist2(h2), v1(p1), v2(p2),axis(a) 
	    {};
   
   ~OlapInfo();
   
   G4bool operator==(const OlapInfo &);
   
   G4NavigationHistory hist1, hist2;
   G4ThreeVector v1, v2;
   G4int axis;
     // geantino ray travelling parallel to that axis detected the olap
   std::vector<OlapStepInfo *> stAB;
   std::vector<OlapStepInfo *> stBA;
   G4String info;
   G4int probNot ;	    
   G4LogicalVolume* originalMother;
};

std::ostream & 
   operator<<(std::ostream& flux,  OlapInfo & oi);


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
