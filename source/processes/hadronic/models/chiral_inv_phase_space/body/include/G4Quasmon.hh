// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Quasmon.hh,v 1.2 1999-12-15 14:52:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4Quasmon_h
#define G4Quasmon_h 1


// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4Quasmon ----------------
//             by Mikhail Kossov, July 1999.
//      class for a Quasmon used by the CHIPS Model
// ------------------------------------------------------------

// Standard G4-headers
#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Parton.hh"
#include "G4PartonVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
//@@ Temporary
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
//CHIPS-headers
//#include "G4QCandidate.hh"
#include "G4QCandidateVector.hh"
#include "G4QResonance.hh"
//#include "G4QResonanceVector.hh"
#include "G4QHadron.hh"
#include "G4QHadronVector.hh"
#include "G4QContent.hh"

class G4Quasmon 
{
public:
  G4Quasmon(G4int Z, G4int N, G4int S,
            G4LorentzVector FourMomentum);                        // Direct Quasmon definition
  G4Quasmon(G4int projPDG, G4int targPDG,                         // Set Valence Quark Content
            G4LorentzVector proj4Mom, G4LorentzVector targ4Mom);  // Interaction definition
  G4Quasmon(const G4Quasmon &right);                              // Quasmon duplication

  ~G4Quasmon();
  int operator==(const G4Quasmon &right) const;
  int operator!=(const G4Quasmon &right) const;

  //@@...Makes sence only for G4NuclearQuasmon (if any)
  //const G4ThreeVector& GetPosition() const;
  //void SetPosition(const G4ThreeVector &aPosition);

  //const G4QHadronVector* GetQHadronList() const;

  G4QHadronVector HadronizeQuasmon();
  G4LorentzVector Get4Momentum() const;

private:  
  G4double   GetQPartonMomentum(G4double kMin);
  void       CalculateContentOfQPartons(G4double qMass);
  void       CalculateHadronizationProbabilities(G4double kVal);
  G4double   GetMinSqThresh (const G4QContent& qCon);
  G4double   GetMidSqThresh (const G4QContent& qCon);
  G4double   GetMaxLSqThresh(const G4QContent& qCon);
  G4double   GetMaxHSqThresh(const G4QContent& qCon);
  G4int      GetNofSqMassSum(const G4QContent& qCon, G4int& sPDG);
  G4bool     Quasmon2HDecay(const G4int& sPDG, G4LorentzVector& r4Mom);
  G4bool     Quasmon2HDecay(const G4int& rPDG, const G4int& sPDG);
  G4bool     Quasmon3HDecay(const G4int& rPDG, const G4int& sPDG, const G4int& tPDG);
private:
  void InitCandidateVector();
  //@@......Make sence only for G4NuclearQuasmon
  //G4ThreeVector   thePosition;  // Position of the Quasmon
  G4QCandidateVector theQCandidates;  // Vector of possible secondary hadrons
  G4QHadronVector    theQResonances;  // Vector of possible residual resonances
  G4QHadronVector    theQHadrons;     // Vector of generated secondary hadrons
  G4LorentzVector    kQParton;        // 4-momentum of QParton kandidate

  G4int           qZ;                 // a#of "protons" in the Quasmon
  G4int           qN;                 // a#of "neutrons" in the Quasmon
  G4int           qS;                 // a#of "lambdas" in the Quasmon
  G4LorentzVector q4Mom;              // 4-momentum of the Quasmon
  G4QContent      valQ;               // valence quark content
  G4int           nOfQ;               // number of quark-partons just to accelerate @@ ??
};

inline int G4Quasmon::operator==(const G4Quasmon &right) const
{
  return this == &right;
}

inline int G4Quasmon::operator!=(const G4Quasmon &right) const
{
  return this != &right;
}

inline G4LorentzVector G4Quasmon::Get4Momentum() const {return q4Mom;}

// Member functions for the output Hadrons

//inline const G4QHadronVector* G4Quasmon::GetQHadronList() const
//{
//  return &theQHadrons;
//}

// @@ Make sence only for G4NuclearQuasmon
//inline const G4ThreeVector& G4Quasmon::GetPosition() const 
//{
//  return thePosition;
//}

// @@ Make sence only for G4NuclearQuasmon
//inline void G4Quasmon::SetPosition(const G4ThreeVector &aPosition)
//{
//  thePosition= aPosition;
//}

#endif



