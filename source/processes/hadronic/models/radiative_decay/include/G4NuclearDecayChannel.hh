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
  
////////////////////////////////////////////////////////////////////////////////
//
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

  G4NuclearDecayChannel (const G4RadioactiveDecayMode &theMode, G4int Verbose) :
    G4GeneralPhaseSpaceDecay (Verbose) {;}
  // default constructor
  //
  G4NuclearDecayChannel (const G4RadioactiveDecayMode &theMode, G4int Verbose,
			 const G4ParticleDefinition *theParentNucleus, 
                         G4double theBR, G4double theQtransition, G4int A,
                         G4int Z, G4double theDaughterExcitation);
  // constructor decay channel with one decay product
  //
  G4NuclearDecayChannel (const G4RadioactiveDecayMode &theMode, G4int Verbose,
			 const G4ParticleDefinition *theParentNucleus,
                         G4double theBR, G4double theQtransition, G4int A,
                         G4int Z, G4double theDaughterExcitation,
			 const G4String theDaughterName1);
  // constructor decay channel with two decay products
  //
  G4NuclearDecayChannel (const G4RadioactiveDecayMode &theMode, G4int Verbose,
                         const G4ParticleDefinition *theParentNucleus,
                         G4double theBR, G4double theFFN,
                         G4double theQtransition, G4int A, G4int Z,
                         G4double theDaughterExcitation,
                         const G4String theDaughterName1,
                         const G4String theDaughterName2);
  // constructor decay channel with three decay product
  //

  ~G4NuclearDecayChannel () {;}
  // destructor
  //
  G4DecayProducts *DecayIt (G4double theParentMass);
  // Returns the decay products
  //
  inline G4RadioactiveDecayMode GetDecayMode () {return decayMode;}
  // Returns the decay mode
  //
  inline G4double GetDaughterExcitation () {return daughterExcitation;}
  // Returns the excitaion energy of the daughter nuclide
  //
  inline G4ParticleDefinition* GetDaughterNucleus () {return daughterNucleus;}
  // Returns the daughter nuclide.
  //
private:
  G4NuclearDecayChannel (const G4String &theName, const G4String &theParentName,
                         G4double theBR, G4int theNumberOfDaughters,
                         const G4String theDaughterName1,
                         const G4String theDaughterName2,
                         const G4String theDaughterName3,
                         const G4String theDaughterName4);

  G4NuclearDecayChannel (const G4String &theParentName,
                         G4double theBR, G4int theNumberOfDaughters,
                         const G4String& theDaughterName1,
                         const G4String& theDaughterName2 = "",
                         const G4String& theDaughterName3 = "");

  G4NuclearDecayChannel (const G4String &theParentName,
                         G4double theParentMass, G4double theBR,
                         G4int theNumberOfDaughters,
                         const G4String& theDaughterName1,
                         const G4String& theDaughterName2 = "",
                         const G4String& theDaughterName3 = "");

  void FillDaughterNucleus (G4int index, G4int A, G4int Z,
                            G4double theDaughterExcitation);

  G4DecayProducts* BetaDecayIt();
  // to replace the ThreeBodyDecayIt() to generate the correct beta spectrum

protected:
  G4RadioactiveDecayMode decayMode;
  static const G4double  pTolerance;
  static const G4double  levelTolerance;
  G4double               daughterExcitation;
  G4int                  daughterA;
  G4int                  daughterZ;
  G4ParticleDefinition  *daughterNucleus;
  G4DynamicParticle     *dynamicDaughter;
  G4double               Qtransition;
  G4double               FermiFN;
  static const G4bool    FermiOn;
};
#endif


