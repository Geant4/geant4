// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QCandidate.hh,v 1.1 1999-11-17 11:04:14 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4QCandidate_h
#define G4QCandidate_h 1

// ----------------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QCandidate ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for Quasmon initiated Candidates used by the CHIPS Model
// ----------------------------------------------------------------------

#include <iostream.h>
//#include <CLHEP/config/CLHEP.h>
#include "globals.hh"
//#include "G4ThreeVector.hh"
#include "G4Nucleus.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleTable.hh"
#include "G4QContent.hh"

class G4QCandidate
{
public:
  G4QCandidate(){PDGencoding=0;}
  //G4QCandidate(G4int PDGcode, G4Nucleus * aNucleus = NULL);
  G4QCandidate(G4int PDGcode);
  G4QCandidate(const G4QCandidate &right);

  ~G4QCandidate();

  const G4QCandidate & operator=(const G4QCandidate &right);

  int operator==(const G4QCandidate &right) const;

  int operator!=(const G4QCandidate &right) const;
      
  G4int GetPDGcode() const;                               // Get PDG code of the Candidate
      
  G4ParticleDefinition * GetDefinition();                 // Get Particle Definition pointer 

  virtual const G4LorentzVector & Get4Momentum() const;   // get 4-mom of the candidate
  virtual void Set4Momentum(const G4LorentzVector& aMom); // set 4-mom of the candidate
      
  G4QContent GetQC();                                     // Get private quark content
  G4double GetMass();                                     // Get a mass of the candidate
  G4double GetMass2();                                    // Get an m^2 value for the candidate
  G4double GetMassChange();                               // @@ Mass change (i.r.t. Initialization)
  void     SetRelProbability(G4double relP);              // Put the relative prob. of hadronization
  G4double GetRelProbability();                           // Get the relative prob. of hadronization
  void     SetIntegProbability(G4double intP);            // Set integrated prob. for randomization
  G4double GetIntegProbability();                         // Get integrated prob. for randomization

  G4int    GetUContent();      
  G4int    GetDContent();      
  G4int    GetSContent();      
  G4int    GetAntiUContent();      
  G4int    GetAntiDContent();      
  G4int    GetAntiSContent();      

  void     SetUContent(G4int nU);      
  void     SetDContent(G4int nD);
  void     SetSContent(G4int nS);
  void     SetAntiUContent(G4int nAU);    
  void     SetAntiDContent(G4int nAD);
  void     SetAntiSContent(G4int nAS);

private:  
  G4int PDGencoding;                    // Code in Particle Data Group
  G4ParticleDefinition * theDefinition; // Pointer to Particle Definition in PDG Table
  G4LorentzVector theMomentum;          // The 4-mom of candidate. In Vacuum = {(0,0,0),M}
  G4double theMass;                     // @@for internal check only
  G4QContent valQ;                      // valence quark content

  G4double relativeProbability;         // relative probability of hadronization
  G4double integralProbability;         //integrated probability for randomization
};

inline int G4QCandidate::operator==(const G4QCandidate &right) const {return this==&right;}	
inline int G4QCandidate::operator!=(const G4QCandidate &right) const {return this!=&right;}

inline G4QContent G4QCandidate::GetQC()                    {return valQ;}
inline G4int    G4QCandidate::GetPDGcode() const           {return PDGencoding;}
inline G4ParticleDefinition *G4QCandidate::GetDefinition() {return theDefinition;}
inline const    G4LorentzVector& G4QCandidate::Get4Momentum()    const {return theMomentum;}
inline void     G4QCandidate::Set4Momentum(const G4LorentzVector& aMom){theMomentum=aMom;}
inline void     G4QCandidate::SetRelProbability(G4double relP)   {relativeProbability=relP;}
inline G4double G4QCandidate::GetRelProbability()          {return relativeProbability;}
inline void     G4QCandidate::SetIntegProbability(G4double intP) {integralProbability=intP;}
inline G4double G4QCandidate::GetIntegProbability()        {return integralProbability;}

inline G4double G4QCandidate::GetMass()          {return theMomentum.m();}
inline G4double G4QCandidate::GetMass2()         {return theMomentum.m2();}
inline G4double G4QCandidate::GetMassChange()    {return theMass - theMomentum.m();}

inline G4int    G4QCandidate::GetUContent()      {return valQ.GetU();}
inline G4int    G4QCandidate::GetDContent()      {return valQ.GetD();}
inline G4int    G4QCandidate::GetSContent()      {return valQ.GetS();}
inline G4int    G4QCandidate::GetAntiUContent()  {return valQ.GetAU();}
inline G4int    G4QCandidate::GetAntiDContent()  {return valQ.GetAD();}
inline G4int    G4QCandidate::GetAntiSContent()  {return valQ.GetAS();}

inline void     G4QCandidate::SetUContent(G4int nU)     {valQ.SetU(nU);}
inline void     G4QCandidate::SetDContent(G4int nD)     {valQ.SetD(nD);}
inline void     G4QCandidate::SetSContent(G4int nS)     {valQ.SetS(nS);}
inline void     G4QCandidate::SetAntiUContent(G4int nAU){valQ.SetAU(nAU);}
inline void     G4QCandidate::SetAntiDContent(G4int nAD){valQ.SetAD(nAD);}
inline void     G4QCandidate::SetAntiSContent(G4int nAS){valQ.SetAS(nAS);}

#endif
