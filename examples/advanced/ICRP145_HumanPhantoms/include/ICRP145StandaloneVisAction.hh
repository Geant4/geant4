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
/// \file visualization/standalone/include/ICRP145StandaloneVisAction.hh
/// \brief Definition of the ICRP145StandaloneVisAction class
//
//

#ifndef ICRP145STANDALONEVISACTION_HH
#define ICRP145STANDALONEVISACTION_HH

#include "G4VUserVisAction.hh"

#include "globals.hh"

class G4UIExecutive;
class TETDetectorConstruction;

class ICRP145StandaloneVisAction: public G4VUserVisAction {
public:
  // sex: male/female: false/true
  ICRP145StandaloneVisAction(G4bool sex, G4UIExecutive* ui);
private:
  void Draw() override;
  TETDetectorConstruction* fTETDetectorConstruction = nullptr;
};

#endif

