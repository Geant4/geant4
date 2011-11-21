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

using namespace std;

G4Allocator<G4IT> aITAllocator;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
///
// Static functions
///
G4IT* GetIT(const G4Track* track)
{
    return (dynamic_cast<G4IT*>(track->GetUserInformation()));
}

G4IT* GetIT(const G4Track& track)
{
    return (dynamic_cast<G4IT*>(track.GetUserInformation()));
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
///
// Constructors / Destructors
///
G4IT::G4IT() : G4VUserTrackInformation("G4IT"),
    fTrack (0),
    fPreviousIT(0), fNextIT(0), fTrackingInformation(new G4TrackingInformation())
{
    fITBox=0;
    fKDNode = 0 ;
    fTrackNode = 0;
    fParentID_A = 0;
    fParentID_B = 0;
}

// Use only by inheriting classes
G4IT::G4IT(const G4IT& /*right*/) : G4VUserTrackInformation("G4IT"),
    fTrack (0),
    fPreviousIT(0), fNextIT(0), fTrackingInformation(new G4TrackingInformation())
{
    fITBox=0;
    fKDNode = 0 ;
    fTrackNode = 0;
    fParentID_A = 0;
    fParentID_B = 0;
}

// Should not be used
G4IT& G4IT::operator=(const G4IT& right)
{
    if(this == &right) return *this;

    fTrack = 0;
    fPreviousIT = 0;
    fNextIT = 0;
    fTrackingInformation = 0;
    fITBox = 0;
    fKDNode = 0 ;
    fTrackNode = 0;
    fParentID_A = 0;
    fParentID_B = 0;
    return *this;
}

G4IT::G4IT(G4Track * aTrack) : G4VUserTrackInformation("G4IT"),
    fPreviousIT(0), fNextIT(0), fTrackingInformation(new G4TrackingInformation())
{
    fITBox = 0;
    fTrack = aTrack;
    fKDNode = 0 ;
    fTrackNode = 0;
    fParentID_A = 0;
    fParentID_B = 0;
    RecordCurrentPositionNTime();
}

void G4IT::TakeOutBox()
{
    if(fITBox)
    {
        fITBox->Extract(this);
    }

    if(fKDNode)
    {
        InactiveNode(fKDNode);
        fKDNode = 0;
    }
}

G4IT::~G4IT()
{
    TakeOutBox();

    if(fTrackingInformation)
    {
        delete fTrackingInformation;
        fTrackingInformation = 0;
    }

    // Note :
    // G4ITTrackingManager will delete fTrackNode.
    // fKDNode will be deleted when the KDTree is rebuilt
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
///
// Methods
///
void G4IT::RecordCurrentPositionNTime()
{
    if(fTrack)
    {
        fTrackingInformation->RecordCurrentPositionNTime(fTrack);
    }
}

G4bool G4IT::operator<(const G4IT& right) const
{
    if(GetITType() == right.GetITType() )
    {
        return  (this->diff(right)) ;
    }
    else
    {
        return (GetITType() < right.GetITType());
    }
    return false;
}

G4bool G4IT::operator==(const G4IT& right) const
{
    if(GetITType() == right.GetITType() )
    {
        return this->equal(right);
    }
    return false;
}

G4bool G4IT::operator!=(const G4IT& right) const
{
    return !(this->operator==(right));
}
