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
        G4ThreeVector Velocity = (1/std::sqrt(BeamMom.mag2() + sqr(KT.GetDefinition()->GetPDGMass())))*BeamMom;
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
