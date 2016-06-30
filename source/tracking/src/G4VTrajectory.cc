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
// $Id: G4VTrajectory.cc 96231 2016-03-30 10:35:46Z gcosmo $
//
// ---------------------------------------------------------------
//
// G4VTrajectory.cc
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@kekvax.kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------

#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4AttCheck.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Colour.hh"

G4VTrajectory::G4VTrajectory() {;}
G4VTrajectory::~G4VTrajectory() {;}

G4bool G4VTrajectory::operator == (const G4VTrajectory& right) const
{
  return (this==&right);
}

void G4VTrajectory::ShowTrajectory(std::ostream& os) const
{
  // Makes use of attribute values implemented in the concrete class.
  // Note: the user needs to follow with new-line or end-of-string,
  // depending on the nature of os.

  std::vector<G4AttValue>* attValues = CreateAttValues();
  const std::map<G4String,G4AttDef>* attDefs = GetAttDefs();

  // Ensure validity...
  if (G4AttCheck(attValues,attDefs).Check("G4VTrajectory::ShowTrajectory")) {
    return;
  }

  os << "Trajectory:";

  std::vector<G4AttValue>::iterator iAttVal;
  for (iAttVal = attValues->begin();
       iAttVal != attValues->end(); ++iAttVal) {
    std::map<G4String,G4AttDef>::const_iterator iAttDef =
      attDefs->find(iAttVal->GetName());
    os << "\n  " << iAttDef->second.GetDesc()
       << " (" << iAttVal->GetName()
       << "): " << iAttVal->GetValue();
  }

  delete attValues;  // AttValues must be deleted after use.

  //Now do trajectory points...
  for (G4int i = 0; i < GetPointEntries(); i++) {

    G4VTrajectoryPoint* aTrajectoryPoint = GetPoint(i);
    attValues = aTrajectoryPoint->CreateAttValues();
    attDefs = aTrajectoryPoint->GetAttDefs();

    // Ensure validity...
    if (G4AttCheck(attValues,attDefs).Check("G4VTrajectory::ShowTrajectory")) {
      return;
    }

    for (iAttVal = attValues->begin();
	 iAttVal != attValues->end(); ++iAttVal) {
      std::map<G4String,G4AttDef>::const_iterator iAttDef =
	attDefs->find(iAttVal->GetName());
      os << "\n    " << iAttDef->second.GetDesc()
	 << " (" << iAttVal->GetName()
	 << "): " << iAttVal->GetValue();
    }

    delete attValues;  // AttValues must be deleted after use.
  }
  os << std::endl;
}

void G4VTrajectory::DrawTrajectory() const
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if (0 != pVVisManager) {
    pVVisManager->DispatchToModel(*this);
  }
}
