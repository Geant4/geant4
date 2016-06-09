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
// $Id: HadrontherapyRunAction.hh,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------

#ifndef HadrontherapyRunAction_h
#define HadrontherapyRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <iostream.h>

// -------------------------------------------------------------
class G4Run;

#ifndef G4NOHIST
namespace AIDA {
 class ITree;
 class IHistogram1D;
} 
#endif

// -------------------------------------------------------------
class HadrontherapyRunAction : public G4UserRunAction
{
  public:
    HadrontherapyRunAction();
   ~HadrontherapyRunAction();
  
G4double energy[50000];
  
public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
G4int GetProva();

  void  EnergyTotSlice(G4int, G4double); 
  void  EnergyTotMarkus();
 
#ifndef G4NOHIST
    AIDA::IHistogram1D* GetHisto(G4int id) {return histo[id];}
#endif

private:
    void bookHisto();

 private:
G4String histName;
    
#ifndef G4NOHIST    
    AIDA::ITree* tree;
    AIDA::IHistogram1D* histo[10];    
#endif
  
G4int NbOfLayer; 
G4int slice;

private:
G4int prova;
G4double energy_dep;
};

#endif

