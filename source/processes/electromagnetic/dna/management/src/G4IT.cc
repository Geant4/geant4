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
// $Id: G4IT.cc 65022 2012-11-12 16:43:12Z gcosmo $
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
    fpTrack (0),
    fpPreviousIT(0), fpNextIT(0),
    fTrackingInformation()
//    fpTrackingInformation(new G4TrackingInformation())
{
    fpITBox=0;
    fpKDNode = 0 ;
    fpTrackNode = 0;
    fParentID_A = 0;
    fParentID_B = 0;
}

// Use only by inheriting classes
G4IT::G4IT(const G4IT& /*right*/) : G4VUserTrackInformation("G4IT"),
    fpTrack (0),
    fpPreviousIT(0), fpNextIT(0),
    fTrackingInformation()
//    fpTrackingInformation(new G4TrackingInformation())
{
    fpITBox=0;
    fpKDNode = 0 ;
    fpTrackNode = 0;
    fParentID_A = 0;
    fParentID_B = 0;
}

// Should not be used
G4IT& G4IT::operator=(const G4IT& right)
{
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "The assignment operator of G4IT should not be used, this feature is not supported."
                         << "If really needed, please contact the developers.";
    G4Exception("G4IT::operator=(const G4IT& right)","G4IT001",FatalException,exceptionDescription);

    if(this == &right) return *this;

    fpTrack = 0;
    fpITBox = 0;
    fpPreviousIT = 0;
    fpNextIT = 0;
    fpKDNode = 0 ;
    fParentID_A = 0;
    fParentID_B = 0;
//    fpTrackingInformation = 0;
    fpTrackNode = 0;

    return *this;
}

G4IT::G4IT(G4Track * aTrack) : G4VUserTrackInformation("G4IT"),
    fpPreviousIT(0), fpNextIT(0),
    fTrackingInformation()
//    fpTrackingInformation(new G4TrackingInformation())
{
    fpITBox = 0;
    fpTrack = aTrack;
    fpKDNode = 0 ;
    fpTrackNode = 0;
    fParentID_A = 0;
    fParentID_B = 0;
    RecordCurrentPositionNTime();
}

void G4IT::TakeOutBox()
{
    if(fpITBox)
    {
        fpITBox->Extract(this);
    }

    if(fpKDNode)
    {
        InactiveNode(fpKDNode);
        fpKDNode = 0;
    }
}

G4IT::~G4IT()
{
    TakeOutBox();

//    if(fpTrackingInformation)
//    {
//        delete fpTrackingInformation;
//        fpTrackingInformation = 0;
//    }

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
    if(fpTrack)
    {
        fTrackingInformation.RecordCurrentPositionNTime(fpTrack);
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
