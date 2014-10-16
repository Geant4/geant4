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

#ifndef G4NuclearDecayChannel_h
#define G4NuclearDecayChannel_h 1

#include "globals.hh"
#include "G4VDecayChannel.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4GeneralPhaseSpaceDecay.hh"
#include "G4RadioactiveDecayMode.hh"
#include "Randomize.hh"


class G4NuclearDecayChannel : public G4GeneralPhaseSpaceDecay
{
  // class description 
  //
  //  G4NuclearDecayChannel is a derived class from G4GeneralPhaseSpaceDecay,
  //  itself a derived class from G4VDecayChannel. It provides the required
  //  decay channels for all nuclear decay modes and through the DecayIt()
  //  member function returns the decay products. 
  // 
  // class description - end
  
  public: // with description

    G4NuclearDecayChannel(const G4RadioactiveDecayMode& theMode, G4int Verbose)
    : G4GeneralPhaseSpaceDecay(Verbose), decayMode(theMode),Qtransition(0) {;}
    // default constructor

    // Decay channel ctor with one decay product
    G4NuclearDecayChannel(const G4RadioactiveDecayMode& theMode, G4int Verbose,
                          const G4ParticleDefinition* theParentNucleus, 
                          const G4double theBR,
                          const G4double theQtransition,
                          const G4int A, const G4int Z, 
                          const G4double theDaughterExcitation);

    G4NuclearDecayChannel(const G4RadioactiveDecayMode& theMode, G4int Verbose,
                          const G4ParticleDefinition* theParentNucleus,
                          G4double theBR, 
                          const G4double theQtransition,
                          const G4int A, const G4int Z,
                          const G4double theDaughterExcitation,
                          const G4String theDaughterName1);
    // constructor decay channel with two decay products
  
    G4NuclearDecayChannel(const G4RadioactiveDecayMode& theMode, G4int Verbose,
                          const G4ParticleDefinition* theParentNucleus,
                          G4double theBR, G4double theFFN,
                          G4bool betaS, G4RandGeneral* randBeta,
                          const G4double theQtransition,
                          const G4int A, const G4int Z,
                          const G4double theDaughterExcitation,
                          const G4String theDaughterName1,
                          const G4String theDaughterName2);
    // constructor decay channel with three decay product

    virtual ~G4NuclearDecayChannel(); 

    // Returns the decay products
    G4DecayProducts* DecayIt(G4double);

    // Set the half-life threshold for isomer production
    void SetHLThreshold(G4double hl) {halflifethreshold = hl;}

    // Enable/disable ICM
    void SetICM(G4bool icm) {applyICM = icm;} 

    // Enable/disable ARM
    void SetARM (G4bool arm) {applyARM = arm;}
 
    inline G4RadioactiveDecayMode GetDecayMode () {return decayMode;}
    // Returns the decay mode

    inline G4double GetDaughterExcitation () {return daughterExcitation;}
    // Returns the excitaion energy of the daughter nuclide

    inline G4ParticleDefinition* GetDaughterNucleus () {return daughterNucleus;}
    // Returns the daughter nuclide.

  private:
    G4NuclearDecayChannel(const G4String& theName, const G4String& theParentName,
                          G4double theBR, G4int theNumberOfDaughters,
                          const G4String theDaughterName1,
                          const G4String theDaughterName2,
                          const G4String theDaughterName3,
                          const G4String theDaughterName4);

    G4NuclearDecayChannel(const G4String& theParentName,
                          G4double theBR, G4int theNumberOfDaughters,
                          const G4String& theDaughterName1,
                          const G4String& theDaughterName2 = "",
                          const G4String& theDaughterName3 = "");

    G4NuclearDecayChannel(const G4String& theParentName,
                          G4double theParentMass, G4double theBR,
                          G4int theNumberOfDaughters,
                          const G4String& theDaughterName1,
                          const G4String& theDaughterName2 = "",
                          const G4String& theDaughterName3 = "");

    void FillDaughterNucleus(G4int index, G4int A, G4int Z,
                             const G4double theDaughterExcitation);

    G4DecayProducts* BetaDecayIt();
    // to replace the ThreeBodyDecayIt() to generate the correct beta spectrum

    // Add copy and assignment prototypes even though no dynamic allocation
    // of class members
    G4NuclearDecayChannel(const G4NuclearDecayChannel &right);
    G4NuclearDecayChannel& operator=(const G4NuclearDecayChannel &right);

  protected:
    //Data members marked with G4ThreadLocal are not invariant among threads
    //Need to make them TLS
    const G4RadioactiveDecayMode decayMode;
    static const G4double pTolerance;
    static const G4double levelTolerance;
    G4double daughterExcitation;
    G4int daughterA;
    G4int daughterZ;
    G4ParticleDefinition* daughterNucleus;  
    static G4ThreadLocal G4DynamicParticle* dynamicDaughter;     
    const G4double Qtransition;
    G4double halflifethreshold;
    G4bool applyICM;
    G4bool applyARM;
    G4RandGeneral* RandomEnergy;   // not dynamically allocated
};
#endif


