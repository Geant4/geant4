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

#ifndef _G4VUSERBROWNIANACTION_HH
#define _G4VUSERBROWNIANACTION_HH
#include "G4Types.hh"
#include "G4Track.hh"
// experimental
class G4VUserBrownianAction
{
 public:
  G4VUserBrownianAction()          = default;
  virtual ~G4VUserBrownianAction() = default;

  //  virtual G4double GetDiffusionCoefficient(G4Material*,
  //                                           G4Molecule*) { return 0;}

  // If returns true: track is killed
  virtual void Transport(G4ThreeVector&, G4Track* pTrack = nullptr) = 0;
  virtual G4double GetDistanceToBoundary(const G4Track&)            = 0;
};

#endif
