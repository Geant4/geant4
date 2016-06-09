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
// $Id: G4FastSimulationManagerProcess.hh,v 1.14 2007/05/11 13:50:20 mverderi Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// 
//---------------------------------------------------------------
//
//  G4FastSimulationManagerProcess.hh
//
//  Description:
//    The process that triggers parameterised simulations  if any.
//
//  History:
//  Feb 98: Parallel geometry sensitivity. MoraDeFreitas.
//  Oct 97: "Fast" replaces "Parameterisation" in class/method names. 
//          (release B.00 for parameterisation). MoraDeFreitas.
//  Aug 97: First implementation. Verderi && MoraDeFreitas.
//  Apr 98: modified for new particle change.  H.Kurashige
//  Oct 06: Move to parallel geometry scheme. M. Verderi
//  Nov 06: name xxx81 is given for this release. "81" will be
//          removed @ next maj. rel. so that this process becomes
//          the standard one.
//  May 07: remove "81" tags, to migrate to 9.0.
//
//---------------------------------------------------------------


#ifndef G4FastSimulationManagerProcess_hh
#define G4FastSimulationManagerProcess_hh

#include "globals.hh"
#include "G4VProcess.hh"
#include "G4FastSimulationManager.hh"
#include "G4Step.hh"
#include "G4Navigator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VParticleChange.hh"


// ---------------------------------------------------------------------
//
//        G4FastSimulationManagerProcess class
//
// ---------------------------------------------------------------------


// Class Description:
// -- G4VProcess providing the interface between the tracking and the fast simulation.
//

class G4FastSimulationManagerProcess : public G4VProcess
{
public:
  
  // -------------------------
  //  Constructor/Destructor:
  // -------------------------
  // -- Constructor for parameterisation in mass geometry
  G4FastSimulationManagerProcess(const G4String&     processName = "G4FastSimulationManagerProcess",
				 G4ProcessType           theType = fParameterisation);
  
  // -- Contructors for parameterisation attached a parallel geometry.
  // -- Can also be used for the mass geometry, providing world volume name.
  // -- World volume specified by name or pointer.
  G4FastSimulationManagerProcess(const G4String&         processName,
				 const G4String&     worldVolumeName,
				 G4ProcessType               theType = fParameterisation);
  G4FastSimulationManagerProcess(const G4String&         processName,
				 G4VPhysicalVolume*      worldVolume,
				 G4ProcessType               theType = fParameterisation);
  
  virtual ~G4FastSimulationManagerProcess();
  
  // -----------------------
  //   User access methods:
  // -----------------------
  G4VPhysicalVolume* GetWorldVolume() const {return fWorldVolume;}
  
  // -- Set new world volume to the process
  void SetWorldVolume(G4String          );
  void SetWorldVolume(G4VPhysicalVolume*);
  
  
  // --------------------------------------------------------------
  //                      Process interface
  // --------------------------------------------------------------
  
  // -- Start/End tracking:
  void StartTracking(G4Track*);
  void   EndTracking();
  
  // -- PostStep methods:
  G4double PostStepGetPhysicalInteractionLength(const G4Track&                track,
						G4double           previousStepSize,
						G4ForceCondition*         condition);
  
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );
  
  
  
  
  
  
  // -- AtRest methods (still there after many years of no use...):
  G4double AtRestGetPhysicalInteractionLength(const G4Track&,
					      G4ForceCondition*);
  
  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&);
  
  
  // -- unused:
  G4double AlongStepGetPhysicalInteractionLength(const G4Track&, G4double, G4double, G4double&, G4GPILSelection*) {return 0.0;}
  G4VParticleChange* AlongStepDoIt(const G4Track& , const G4Step&) {return 0;}
  
  
  // -- debug:
  void Verbose() const;
  
  
private:
  //-- would be better to my taste to have "const G4VPhysicalVolume* fWorldVolume;", but clashes at compilation
  G4VPhysicalVolume*                 fWorldVolume;
  
  G4bool                          fIsTrackingTime;
  G4bool                             fIsFirstStep;
  G4Navigator*                         fNavigator;
  G4int                           fNavigatorIndex;
  
  
  G4FastSimulationManager* fFastSimulationManager;
  G4bool                   fFastSimulationTrigger;
};

#endif
