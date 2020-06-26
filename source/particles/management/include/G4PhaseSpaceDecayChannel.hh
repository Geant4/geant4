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
// G4PhaseSpaceDecayChannel
//
// Class description:
//
// Concrete decay channel for a particle.

// Author: H.Kurashige, 27 July 1996
// --------------------------------------------------------------------
#ifndef G4PhaseSpaceDecayChannel_hh
#define G4PhaseSpaceDecayChannel_hh 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"
#include "G4Cache.hh"

class G4PhaseSpaceDecayChannel : public G4VDecayChannel
{
  public:

    enum { MAX_N_DAUGHTERS=4 }; 

    G4PhaseSpaceDecayChannel(G4int Verbose = 1);
    G4PhaseSpaceDecayChannel(const G4String& theParentName,
                                   G4double  theBR,
                                   G4int     theNumberOfDaughters,
                             const G4String& theDaughterName1,
                             const G4String& theDaughterName2 = "",
                             const G4String& theDaughterName3 = "",
                             const G4String& theDaughterName4 = "" );
      // Constructors 

    virtual ~G4PhaseSpaceDecayChannel();
      // Destructor

    G4bool SetDaughterMasses( G4double masses[]);
      // Give daughter masses instead of sampling masses 
      // according to PDG mass and width
  
    G4bool SampleDaughterMasses();

    virtual G4DecayProducts* DecayIt(G4double);   
    G4bool IsOKWithParentMass(G4double parentMass);

    static G4double Pmx(G4double e, G4double p1, G4double p2);
 
  private:

    G4DecayProducts* OneBodyDecayIt();
    G4DecayProducts* TwoBodyDecayIt();
    G4DecayProducts* ThreeBodyDecayIt();
    G4DecayProducts* ManyBodyDecayIt();
     
    G4Cache<G4double> current_parent_mass;  // A thread-local object
    G4double givenDaughterMasses[MAX_N_DAUGHTERS];
    G4bool useGivenDaughterMass = false;
};  

#endif
