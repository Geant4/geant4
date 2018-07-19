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
// $Id: G4PhaseSpaceDecayChannel.hh 95906 2016-03-02 10:56:50Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      27 July 1996 H.Kurashige
//      30 May 1997 H.Kurashige
//      24 May 2012 H.Kurashige    Add current_parent_mass;
// ------------------------------------------------------------
#ifndef G4PhaseSpaceDecayChannel_h
#define G4PhaseSpaceDecayChannel_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"
#include "G4Cache.hh"

class G4PhaseSpaceDecayChannel :public G4VDecayChannel
{
  public:  // With Description
    enum {MAX_N_DAUGHTERS=4}; 
    //Constructors 
      G4PhaseSpaceDecayChannel(G4int Verbose = 1);
      G4PhaseSpaceDecayChannel(const G4String& theParentName,
			       G4double        theBR,
			       G4int           theNumberOfDaughters,
			       const G4String& theDaughterName1,
			       const G4String& theDaughterName2 = "",
			       const G4String& theDaughterName3 = "",
			       const G4String& theDaughterName4 = ""   );

   // give daughter masses instead of sampling masses 
   // according to PDG mass and width. 
     G4bool SetDaughterMasses( G4double masses[]);
  
     G4bool SampleDaughterMasses();

  public: 
   //  Destructor
     virtual ~G4PhaseSpaceDecayChannel();

  public:  // With Description
     virtual G4DecayProducts *DecayIt(G4double);   
     G4bool  IsOKWithParentMass(G4double parentMass);

  public: 
     static G4double Pmx(G4double e, G4double p1, G4double p2);
 
  private:
     //A thread-local object
     G4Cache<G4double> current_parent_mass;
     G4DecayProducts *OneBodyDecayIt();
     G4DecayProducts *TwoBodyDecayIt();
     G4DecayProducts *ThreeBodyDecayIt();
     G4DecayProducts *ManyBodyDecayIt();
     G4bool           useGivenDaughterMass;
     G4double         givenDaughterMasses[MAX_N_DAUGHTERS];
};  


#endif
