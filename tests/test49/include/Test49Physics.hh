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
// -------------------------------------------------------------
// Short description of the G4 Class made by M.K. (asFarAsHeNnderstands)
// This is a substitute for the Physics List in Test49
// =====================================================================
//
//      ---------- Test49Physics -------
// 
//    Converted from Test29 to Test49 by Mikhail Kossov, 29 Jan 2005 
//
//---------------------------------------------------------------------------------------

#ifndef Test49Physics_h
#define Test49Physics_h 1

#include "globals.hh"
#include "Test49HadronProduction.hh"
#include "G4VRestProcess.hh"
#include "G4VDiscreteProcess.hh"
#include "G4QCaptureAtRest.hh"
#include "G4StringChipsInterface.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4BinaryCascade.hh"

class G4VProcess;
class G4Material;

class Test49Physics
{
public:
  Test49Physics();
  ~Test49Physics();

  G4VProcess* GetProcess(const G4String&, const G4String&, G4Material*, G4bool, G4bool);
  G4double GetNucleusMass() {return theHDProcess->GetMass();}; //Only for HadronicProcesses
  G4ExcitationHandler* GetDeExcitation() {return theDeExcitation;};//Only for HadrProcesses
  G4PreCompoundModel* GetPreCompound() {return thePreCompound;};// Only for HadronProcesses
  //void setCutOnP(G4double val) {if(hkmod) hkmod->setCutOnP(val);};
  //void setCutOnPPP(G4double val) {if(hkmod) hkmod->setCutOnPPP(val);};

private:
  // Body
  G4VDiscreteProcess*     theWDProcess; // ElectroWeAk(W) Discrete(D) Process (onFlight)
  G4VRestProcess*         theWRProcess; // ElectroWeak(W) Rest(R) Process (atRest)
  Test49HadronProduction* theHDProcess; // Hadronic(H) Discrete(D) Process (onFlight)
  G4ExcitationHandler*    theDeExcitation;
  G4PreCompoundModel*     thePreCompound;
  G4BinaryCascade*        hkmod;
};

#endif

 


