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
// $Id: ParExN04MarshaledHits.cc,v 1.1 2002-03-05 15:22:18 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
//                   Parallel Library for Geant4
//
//             Gene Cooperman <gene@ccs.neu.edu>, 2001
// --------------------------------------------------------------------
// In this version, we MarshalHit and UnmarshalHit for
// TMarshaledHitsCollection<ExN04CalorimeterHit>
// more efficiently by ignoring position, rotation matrix,
// and pLogV under assumption that they won't be needed for analysis.
// =====================================================================
// Must define MarshalHit and UnmarshalHit for each derived class of G4VHit
// MarshalHit() and UnmarshalHit() have defaults given by the member
// function of TMarshaledHitsCollection<T>.
// The default copies the entire class into a single buffer.  If your
// hit class includes pointers, and if you use those pointer in your
// analysis (AnalyzeEvent), then you will need to override the default
// with an explicit specialization that marshals the pointers.
// You may also wish to override the default in order to write more
// efficient marshaling code.
// In this case, we choose to marshal only the essential fields of
// ExN04CalorimeterHit.
// --------------------------------------------------------------------

#include "ParExMarshaledHits.hh"
#include "ParMarshaledObj.hh"
#include "ExN04TrackerHit.hh"
#include "ExN04MuonHit.hh"
#include "ExN04CalorimeterHit.hh"

template <>
inline void TMarshaledHitsCollection<ExN04CalorimeterHit>::MarshalHit( ExN04CalorimeterHit *aHit )
{ 
  MarshaledObj::Marshal( aHit->GetZ() );
  MarshaledObj::Marshal( aHit->GetPhi() );
  MarshaledObj::Marshal( aHit->GetEdep() );
}
template <>
inline void TMarshaledHitsCollection<ExN04CalorimeterHit>::UnmarshalHit( ExN04CalorimeterHit *aHit )
{
  if (isNewHitsCollection) {
    G4int ZCellID = UnmarshalAsG4int();
    G4int PhiCellID = UnmarshalAsG4int();
    aHit->SetCellID( ZCellID, PhiCellID );
    aHit->SetEdep( UnmarshalAsG4double() );
    // Here, we decide that G4ThreeVector and G4RotationMatrix and pLogV
    //   will not be needed for analysis.
  }
  else
    G4Exception(
      "Event Level Parallelism:  Every hit collection sent to master should be"
      " new.\n  Current hit collection, " + collectionName + ", is not new.\n"
    );
}

// =====================================================================
// Dispatch marshaling and unmarshaling of HC's according to type of G4VHit
//   When unmarshaling, we will see the "G4String HCname", and have
//   to create a HC of the correct type.
//   So, You must include an "if" or "else if" for each type of hit.

// NOTE TO MYSELF:  Is there a cleaner way to do this?
//   If types were objects, we'd just set up a table here:
//     { {"trackerCollection", TMarshaledHitsCollection<ExN02TrackerHit>},
//       ... }
//   and then write general code to use the table.

void MarshaledHCofThisEvent::MarshalHitsCollection
				(G4VHitsCollection * aBaseHC, G4String HCname)
{ if (HCname == "trackerCollection")
    ((TMarshaledHitsCollection<ExN04TrackerHit> *)this)
      -> Marshal( aBaseHC );
  else if (HCname == "muonCollection")
    ((TMarshaledHitsCollection<ExN04MuonHit> *)this)
      -> Marshal( aBaseHC );
  else if (HCname == "calCollection")
    ((TMarshaledHitsCollection<ExN04CalorimeterHit> *)this)
      -> Marshal( aBaseHC );
  else G4Exception( "TMarshaledHitsCollection<T>::MarshaledHCofThisEvent :"
		    " HCname, " + HCname + ", not found." );
}

void MarshaledHCofThisEvent::UnmarshalSlaveHitsCollection(G4String & HCname, G4String & SDname)
{ if (HCname == "trackerCollection")
    ((TMarshaledHitsCollection<ExN04TrackerHit> *)this)
      -> UnmarshalSlaveHitsCollection( HCname, SDname );
  else if (HCname == "muonCollection")
    ((TMarshaledHitsCollection<ExN04MuonHit> *)this)
      -> UnmarshalSlaveHitsCollection( HCname, SDname );
  else if (HCname == "calCollection")
    ((TMarshaledHitsCollection<ExN04CalorimeterHit> *)this)
      -> UnmarshalSlaveHitsCollection( HCname, SDname );
  else G4Exception( "TMarshaledHitsCollection<T>::UnmarshalSlaveHitsCollection"
                    " : HCname, " + HCname + ", not found." );
}
