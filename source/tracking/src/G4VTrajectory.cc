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
// $Id: G4VTrajectory.cc,v 1.8 2005/11/15 03:52:26 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
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
    os << "\n  " << iAttDef->second.GetDesc() << ": "
       << iAttVal->GetValue();
  }

  delete attValues;  // AttValues must be deleted after use.

  //Now do trajectory points...
  for (G4int i = 0; i < GetPointEntries(); i++) {

    G4VTrajectoryPoint* aTrajectoryPoint = GetPoint(i);
    std::vector<G4AttValue>* attValues
      = aTrajectoryPoint->CreateAttValues();
    const std::map<G4String,G4AttDef>* attDefs
      = aTrajectoryPoint->GetAttDefs();

    // Ensure validity...
    if (G4AttCheck(attValues,attDefs).Check("G4VTrajectory::ShowTrajectory")) {
      return;
    }

    std::vector<G4AttValue>::iterator iAttVal;
    for (iAttVal = attValues->begin();
	 iAttVal != attValues->end(); ++iAttVal) {
      std::map<G4String,G4AttDef>::const_iterator iAttDef =
	attDefs->find(iAttVal->GetName());
      os << "\n    " << iAttDef->second.GetDesc() << ": "
	 << iAttVal->GetValue();
    }

    delete attValues;  // AttValues must be deleted after use.
  }
}

void G4VTrajectory::DrawTrajectory(G4int i_mode) const
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if (0 != pVVisManager) {
    pVVisManager->DispatchToModel(*this, i_mode);
  }
}
