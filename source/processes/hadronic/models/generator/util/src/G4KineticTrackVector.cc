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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "G4KineticTrackVector.hh"

G4KineticTrackVector::G4KineticTrackVector()
{
}


//****************************************************************************************************************
// These methods were implemented by Maxim Komogorov
// Maxim.Komogorov@cern.ch

void G4KineticTrackVector::BoostBeam(G4ThreeVector& BeamMom)
    {
    for(unsigned int c1 = 0; c1 < size(); c1++)
        {
        G4KineticTrack& KT =**(begin()+c1);
        G4LorentzVector Mom = KT.Get4Momentum();        
        G4ThreeVector Velocity = (1/sqrt(BeamMom.mag2() + sqr(KT.GetDefinition()->GetPDGMass())))*BeamMom;
        Mom.boost(Velocity);
        KT.Set4Momentum(Mom);
        }
    }
//--------------------------------------------------------------------------------------------------------------

void G4KineticTrackVector::Boost(G4ThreeVector& Velocity)
    { 
    for(unsigned int c1 = 0; c1 < size(); c1++)
        {
        G4KineticTrack& KT =**(begin()+c1);
        G4LorentzVector Mom = KT.Get4Momentum();        
        Mom.boost(Velocity);
        KT.Set4Momentum(Mom);
        }
    }

//--------------------------------------------------------------------------------------------------------------

void G4KineticTrackVector::Shift(G4ThreeVector& Pos)
    {
    for(unsigned int c1 = 0; c1 < size(); c1++)
        {
        G4KineticTrack& KT =**(begin()+c1);
        KT.SetPosition(KT.GetPosition() + Pos);
        }
    }
 

//****************************************************************************************************************
