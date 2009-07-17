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
// $Id: HadrontherapyDetectorMessenger.hh;
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
#ifndef IAEADetectorMessenger_h
#define IAEADetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class IAEADetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

class IAEADetectorMessenger: public G4UImessenger
{
  public:
  IAEADetectorMessenger(IAEADetectorConstruction* );
  ~IAEADetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
private:

  // Pointer to the detector component
  IAEADetectorConstruction* IAEADetector;

  G4UIcmdWithADoubleAndUnit* setIAEAWaterPhantomThicknessCmd;
  //< UI command to set the thickness of the water phantom (also moves other objects accordingly)
};
#endif

