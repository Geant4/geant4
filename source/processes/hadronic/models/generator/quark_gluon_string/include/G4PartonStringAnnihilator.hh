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
// $Id: G4PartonStringAnnihilator.hh,v 1.6 2002/12/12 19:17:34 gunter Exp $
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 9-Oct-1998
// -----------------------------------------------------------------------------
#ifndef G4PartonStringAnnihilator_h
#define G4PartonStringAnnihilator_h 1

#include "G4KineticTrack.hh"
#include "G4ExcitedString.hh"
#include "G4MesonSplitter.hh"
#include "G4BaryonSplitter.hh"


class G4PartonStringAnnihilator
    {   
public:
    G4PartonStringAnnihilator();
    
public:
   
    G4double GetCrossSection(G4KineticTrack& aTarget, G4KineticTrack& aProjectile);
    G4ExcitedString* GetString(G4KineticTrack& Target, G4KineticTrack& Projectile);

    G4bool SplitUpBarion(G4int Encoding, G4int* q_or_qqbar, G4int* qbar_or_qq); //!
 
    G4bool isMeson(G4int Encoding); 
    G4bool isBarion(G4int Encoding);
    G4bool isAntiParticle(G4int Encoding);

private:    
    G4bool GetStringEnds(G4int Projectile, G4int Target, G4int* Left, G4int* Right);
    G4bool FindDiquark(G4int Encoding, G4int Quark, G4int* Diquark);
    G4bool FindQuark(G4int Encoding, G4int Diquark, G4int* Quark);

private:

    G4double widthOfPtSquare;
    G4MesonSplitter theMesonSplitter;
    G4BaryonSplitter theBaryonSplitter;
    };

// *************************************************************************************************
#endif
