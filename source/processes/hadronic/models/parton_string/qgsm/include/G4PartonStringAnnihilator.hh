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
// $Id: G4PartonStringAnnihilator.hh,v 1.2 2006/06/29 20:56:11 gunter Exp $
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
