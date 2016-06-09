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
// $Id: MedLinacEventAction.hh,v 1.4 2005/07/03 23:27:36 mpiergen Exp $
//
//
// Code developed by: M. Piergentili
//
 
#ifndef MedLinacEventAction_h
#define MedLinacEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

//class G4Event;
//class MedLinacDetectorConstruction; 
//class MedLinacAnalysisManager;

//*********************************************************************

class MedLinacEventAction : public G4UserEventAction
{
  public:
    MedLinacEventAction();
   ~MedLinacEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

private:
  G4String drawFlag; //Visualisation flag
  };

//*********************************************************************

#endif

    
