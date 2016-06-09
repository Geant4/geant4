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
// $Id: HadrontherapyRunAction.hh,v 3.0, September 2004;
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#ifndef HadrontherapyRunAction_h
#define HadrontherapyRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4RunManager.hh"
#include "globals.hh"

class G4Run;
class HadrontherapyAnalysisManager;
class HadrontherapyDetectorConstruction;
class HadrontherapyRunMessenger;
class HadrontherapyFactory;
class HadrontherapyFactoryIr;
class HadrontherapyFactoryI;



class HadrontherapyRunAction : public G4UserRunAction
{
public:
  HadrontherapyRunAction();
  ~HadrontherapyRunAction();

public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run* );
  void SelectEnergy(G4int); 

private:
  
  G4String sensitiveDetectorName;
  HadrontherapyDetectorConstruction* detector;
  HadrontherapyRunMessenger* runMessenger;
  HadrontherapyFactory *factory; 
  G4int sourceChoice; //select primary particle 

 
  G4double matrix[80][80][80];
  G4int indexI; 
  G4int indexJ; 
  G4int indexK;
};
#endif



