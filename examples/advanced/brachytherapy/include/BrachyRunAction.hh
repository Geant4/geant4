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
// $Id: BrachyRunAction.hh,v 1.6 2002-11-18 15:18:37 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//
//    ********************************************
//    *                                          *
//    *      BrachyRunAction.hh                  *
//    *                                          *
//    ********************************************
//
#ifndef BrachyRunAction_h
#define BrachyRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4RunManager.hh"
#include "globals.hh"


class G4Run;
class BrachyAnalysisManager;
class BrachyDetectorConstruction;
class BrachyRunMessenger;

class BrachyFactory;
class BrachyFactoryIr;
class BrachyFactoryI;
class BrachyRunAction : public G4UserRunAction
{
  public:
    BrachyRunAction(G4String& );
   ~BrachyRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run* );
  void SelectEnergy(G4int); 
 private:
  
   G4String      SDname;
   BrachyFactory  *factory; 
   G4VUserDetectorConstruction* pDetector;
   BrachyDetectorConstruction* pDet;
 
  G4int a;
  BrachyRunMessenger* pRun;
 
};

#endif



