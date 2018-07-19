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
// $Id: G4IT.cc 103042 2017-03-10 11:50:07Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4IT.hh"
#include "G4KDTree.hh"
#include "G4ITBox.hh"
#include "G4Track.hh"
#include "G4TrackList.hh"
#include "G4TrackingInformation.hh"

using namespace std;

//------------------------------------------------------------------------------
//
// Static functions
//
G4IT* GetIT(const G4Track* track)
{
  return (dynamic_cast<G4IT*>(track->GetUserInformation()));
}

G4IT* GetIT(const G4Track& track)
{
  return (dynamic_cast<G4IT*>(track.GetUserInformation()));
}

template<>
G4KDNode<G4IT>::~G4KDNode(){
  fPoint->SetNode(nullptr);
}

//------------------------------------------------------------------------------
//
// Constructors / Destructors
//
G4IT::G4IT() :
    G4VUserTrackInformation("G4IT"),
    fpTrack(nullptr),
    fpPreviousIT(nullptr),
    fpNextIT(nullptr),
    fpTrackingInformation(new G4TrackingInformation())
{
  fpITBox = nullptr;
  fpKDNode = nullptr;
  fpTrackNode = nullptr;
  fParentID_A = 0;
  fParentID_B = 0;
}

// Use only by inheriting classes
G4IT::G4IT(const G4IT& /*right*/) :
    G4VUserTrackInformation("G4IT"),
    fpTrack(nullptr),
    fpPreviousIT(nullptr),
    fpNextIT(nullptr),
    fpTrackingInformation(new G4TrackingInformation())
{
  fpITBox = nullptr;
  fpKDNode = nullptr;
  fpTrackNode = nullptr;
  fParentID_A = 0;
  fParentID_B = 0;
}

// Should not be used
G4IT& G4IT::operator=(const G4IT& right)
{
  G4ExceptionDescription exceptionDescription;
  exceptionDescription
      << "The assignment operator of G4IT should not be used, "
          "this feature is not supported."
      << "If really needed, please contact the developers.";
  G4Exception("G4IT::operator=(const G4IT& right)",
              "G4IT001",
              FatalException,
              exceptionDescription);

  if (this == &right) return *this;

  fpTrack = nullptr;
  fpITBox = nullptr;
  fpPreviousIT = nullptr;
  fpNextIT = nullptr;
  fpKDNode = nullptr;
  fParentID_A = 0;
  fParentID_B = 0;
  fpTrackingInformation = nullptr;
  fpTrackNode = nullptr;

  return *this;
}

G4IT::G4IT(G4Track * aTrack) :
    G4VUserTrackInformation("G4IT"),
    fpPreviousIT(0),
    fpNextIT(0),
    fpTrackingInformation(new G4TrackingInformation())
{
  fpITBox = 0;
  fpTrack = aTrack;
  fpKDNode = nullptr;
  fpTrackNode = nullptr;
  fParentID_A = 0;
  fParentID_B = 0;
  RecordCurrentPositionNTime();
}

void G4IT::TakeOutBox()
{
  if(fpITBox)
  {
    fpITBox->Extract(this);
    fpITBox = nullptr;
  }

  if(fpTrackNode)
  {
    delete fpTrackNode;
    fpTrackNode = nullptr;
  }

  if(fpKDNode)
  {
    InactiveNode(fpKDNode);
    fpKDNode = nullptr;
  }
}

G4IT::~G4IT()
{
  TakeOutBox();

  if(fpTrackingInformation)
  {
    delete fpTrackingInformation;
    fpTrackingInformation = nullptr;
  }

// Note :
// G4ITTrackingManager will delete fTrackNode.
// fKDNode will be deleted when the KDTree is rebuilt
}

//------------------------------------------------------------------------------
///
// Methods
///

G4bool G4IT::operator<(const G4IT& right) const
{
  if (GetITType() == right.GetITType())
  {
    return (this->diff(right));
  }
  else
  {
    return (GetITType() < right.GetITType());
  }
  return false;
}

G4bool G4IT::operator==(const G4IT& right) const
{
  if (GetITType() == right.GetITType())
  {
    return this->equal(right);
  }
  return false;
}

G4bool G4IT::operator!=(const G4IT& right) const
{
  return !(this->operator==(right));
}

double G4IT::operator[](int i) const
{
  return fpTrack->GetPosition()[i];
}

//------------------------------------------------------------------------------

const G4ThreeVector& G4IT::GetPosition() const
{
  if (fpTrack) return GetTrack()->GetPosition();
  return *(new G4ThreeVector());
}

void G4IT::RecordCurrentPositionNTime()
{
  if (fpTrack)
  {
    fpTrackingInformation->RecordCurrentPositionNTime(fpTrack);
  }
}

G4double G4IT::GetPreStepGlobalTime() const
{
  return fpTrackingInformation->GetPreStepGlobalTime();
}

G4double G4IT::GetPreStepLocalTime() const
{
  return fpTrackingInformation->GetPreStepLocalTime();
}

const G4ThreeVector& G4IT::GetPreStepPosition() const
{
  return fpTrackingInformation->GetPreStepPosition();
}

