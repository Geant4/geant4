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
// $Id: G4VTrajectory.cc,v 1.2 2002-12-11 15:45:03 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Colour.hh"

void G4VTrajectory::ShowTrajectory(G4std::ostream& os) const
{
  // Makes use of attribute values implemented in the concrete class.
  // Note: the user needs to follow with new-line or end-of-string,
  // depending on the nature of os.

  G4std::vector<G4AttValue>* attValues = CreateAttValues();

  if (!attValues) {
    os << "G4VTrajectory::ShowTrajectory: no attribute values defined.";
    return;
  }

  const G4std::map<G4String,G4AttDef>* attDefs = GetAttDefs();
  if (!attDefs) {
    os << "G4VTrajectory::ShowTrajectory:"
      "\n  ERROR: no attribute definitions for attribute values.";
    return;
  }

  os << "Trajectory:";

  G4std::vector<G4AttValue>::iterator iAttVal;
  for (iAttVal = attValues->begin();
       iAttVal != attValues->end(); ++iAttVal) {
    G4std::map<G4String,G4AttDef>::const_iterator iAttDef =
      attDefs->find(iAttVal->GetName());
    if (iAttDef == attDefs->end()) {
      os << "G4VTrajectory::ShowTrajectory:"
	"\n  WARNING: no matching definition for attribute \""
	 << iAttVal->GetName() << "\", value: "
	 << iAttVal->GetValue();
    }
    else {
      os << "\n  " << iAttDef->second.GetDesc() << ": "
	 << iAttVal->GetValue();
    }
  }

  delete attValues;  // AttValues must be deleted after use.

  //Now do trajectory points...
  for (G4int i = 0; i < GetPointEntries(); i++) {
    G4VTrajectoryPoint* aTrajectoryPoint = GetPoint(i);
    G4std::vector<G4AttValue>* attValues
      = aTrajectoryPoint->CreateAttValues();
    if (!attValues) {
      os <<
	"\nG4VTrajectory::ShowTrajectory: no attribute values"
	" for trajectory point defined.";
    }
    else {
      const G4std::map<G4String,G4AttDef>* attDefs
	= aTrajectoryPoint->GetAttDefs();
      if (!attDefs) {
	os << "\nG4VTrajectory::ShowTrajectory:"
	  "\n  ERROR: no attribute definitions for attribute values"
	  " for trajectory point defined.";
      }
      else {
	G4std::vector<G4AttValue>::iterator iAttVal;
	for (iAttVal = attValues->begin();
	     iAttVal != attValues->end(); ++iAttVal) {
	  G4std::map<G4String,G4AttDef>::const_iterator iAttDef =
	    attDefs->find(iAttVal->GetName());
	  if (iAttDef == attDefs->end()) {
	    os << "\nG4VTrajectory::ShowTrajectory:"
	      "\n  WARNING: no matching definition for trajectory"
	      " point attribute \""
	       << iAttVal->GetName() << "\", value: "
	       << iAttVal->GetValue();
	  }
	  else {
	    os << "\n    " << iAttDef->second.GetDesc() << ": "
	       << iAttVal->GetValue();
	  }
	}
      }
      delete attValues;  // AttValues must be deleted after use.
    }
  }
}

void G4VTrajectory::DrawTrajectory(G4int i_mode) const
{
  // If i_mode>=0, draws a trajectory as a polyline (blue for
  // positive, red for negative, green for neutral) and, if i_mode!=0,
  // adds markers - yellow circles for step points and magenta squares
  // for auxiliary points, if any - whose screen size in pixels is
  // given by abs(i_mode)/1000.  E.g: i_mode = 5000 gives easily
  // visible markers.

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (!pVVisManager) return;

  const G4double markerSize = G4std::abs(i_mode)/1000;
  G4bool lineRequired (i_mode >= 0);
  G4bool markersRequired (markerSize > 0.);

  G4Polyline trajectoryLine;
  G4Polymarker stepPoints;
  G4Polymarker auxiliaryPoints;

  for (G4int i = 0; i < GetPointEntries() ; i++) {
    G4VTrajectoryPoint* aTrajectoryPoint = GetPoint(i);
    const G4std::vector<G4ThreeVector>* auxiliaries
      = aTrajectoryPoint->GetAuxiliaryPoints();
    if (auxiliaries) {
      for (size_t iAux = 0; iAux < auxiliaries->size(); ++iAux) {
	const G4ThreeVector pos((*auxiliaries)[iAux]);
	if (lineRequired) {
	  trajectoryLine.push_back(pos);
	}
	if (markersRequired) {
	  auxiliaryPoints.push_back(pos);
	}
      }
    }
    const G4ThreeVector pos(aTrajectoryPoint->GetPosition());
    if (lineRequired) {
      trajectoryLine.push_back(pos);
    }
    if (markersRequired) {
      stepPoints.push_back(pos);
    }
  }

  if (lineRequired) {
    G4Colour colour;
    const G4double charge = GetCharge();
    if(charge>0.)      colour = G4Colour(0.,0.,1.); // Blue = positive.
    else if(charge<0.) colour = G4Colour(1.,0.,0.); // Red = negative.
    else               colour = G4Colour(0.,1.,0.); // Green = neutral.
    G4VisAttributes trajectoryLineAttribs(colour);
    trajectoryLine.SetVisAttributes(&trajectoryLineAttribs);
    pVVisManager->Draw(trajectoryLine);
  }
  if (markersRequired) {
    auxiliaryPoints.SetMarkerType(G4Polymarker::squares);
    auxiliaryPoints.SetScreenSize(markerSize);
    auxiliaryPoints.SetFillStyle(G4VMarker::filled);
    G4VisAttributes auxiliaryPointsAttribs(G4Colour(0.,1.,1.));  // Magenta
    auxiliaryPoints.SetVisAttributes(&auxiliaryPointsAttribs);
    pVVisManager->Draw(auxiliaryPoints);

    stepPoints.SetMarkerType(G4Polymarker::circles);
    stepPoints.SetScreenSize(markerSize);
    stepPoints.SetFillStyle(G4VMarker::filled);
    G4VisAttributes stepPointsAttribs(G4Colour(1.,1.,0.));  // Yellow.
    stepPoints.SetVisAttributes(&stepPointsAttribs);
    pVVisManager->Draw(stepPoints);
  }
}
