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
// $Id: ParExN02MarshaledHits.cc,v 1.1 2002-03-05 15:22:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
//                   Parallel Library for Geant4
//
//             Gene Cooperman <gene@ccs.neu.edu>, 2001
// --------------------------------------------------------------------
// Customize this file for your application.
// --------------------------------------------------------------------

#include "ParExMarshaledHits.hh"
#include "ParMarshaledObj.hh"
#include "ExN02TrackerHit.hh"

// =====================================================================
// Must define MarshalHit and UnmarshalHit for each derived class of G4VHit

// MarshalHit() and UnmarshalHit() have defaults given by the member
//   function of TMarshaledHitsCollection<T>.
// The default copies the entire class into a single buffer.  If your
//   hit class includes pointers, and if you use those pointer in your
//   analysis (AnalyzeEvent), then you will need to override the default
//   with an explicit specialization that marshals the pointers.
//   You may also wish to override the default in order to write
//   more efficient marshaling code.

// See ParExMarshaledHits.cc-N04-efficient for an example of overriding
//   the default.

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
    ((TMarshaledHitsCollection<ExN02TrackerHit> *)this)
      -> Marshal( aBaseHC );
  else G4Exception( "TMarshaledHitsCollection<T>::MarshaledHCofThisEvent :"
		    " HCname, " + HCname + ", not found." );
}

void MarshaledHCofThisEvent::UnmarshalSlaveHitsCollection(G4String & HCname, G4String & SDname)
{ if (HCname == "trackerCollection")
    ((TMarshaledHitsCollection<ExN02TrackerHit> *)this)
       -> UnmarshalSlaveHitsCollection( HCname, SDname );
  else G4Exception( "TMarshaledHitsCollection<T>::UnmarshalSlaveHitsCollection"
                    " : HCname, " + HCname + ", not found." );
}
