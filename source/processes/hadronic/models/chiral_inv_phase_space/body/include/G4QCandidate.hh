// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QCandidate.hh,v 1.3 2000-08-16 07:32:45 hpw Exp $
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

#include "G4QNucleus.hh"
#include "G4QParentClusterVector.hh"

class G4QCandidate : public G4QHadron
{
public:
  G4QCandidate();                                            // Default Constructor
  G4QCandidate(G4int PDGcode);                               //Constructor
  G4QCandidate(const G4QCandidate &right);                   // Copy Constructor
  ~G4QCandidate();                                           // Destructor
  // Overloaded Operators
  const G4QCandidate & operator=(const G4QCandidate &right);
  G4int operator==(const G4QCandidate &right) const;
  G4int operator!=(const G4QCandidate &right) const;
  // Specific Selectors
  G4QParentCluster* TakeParClust(G4int nPC);// Get pointer to Parent Cluster from ParClastVector
  G4int    GetPClustEntries()    const;     // Get a number of Parent Clusters in ParClastVector
  G4bool   GetPossibility()      const;     // Get possibility(true)/forbiding(false) to be a hadr./fragm.
  G4bool   GetParPossibility()   const;     // Get possibility(true)/forbiding(false) to be a parent
  G4double GetKMin()             const;     // Get k-minimal for the candidate
  G4double GetDenseProbability() const;     // Get dense-probability for the candidate 
  G4double GetPreProbability()   const;     // Get pre-probability for the candidate 
  G4double GetRelProbability()   const;     // Get the relative probility of hadronization
  G4double GetIntegProbability() const;     // Get integrated probability for randomization
  G4double GetSecondRelProb()    const;     // Get 2nd relative probility of hadronization
  G4double GetSecondIntProb()    const;     // Get 2nd integ. probability for randomization
  // Specific Modifiers
  void ClearParClustVector();               // Clear the Parent Claster Vector of the Candidate
  void FillPClustVec(G4QParentCluster* pCl);// Fill pointer to Parent Cluster in ParClastVector
  void SetPossibility(G4bool choice);       // Set possibility(true)/forbiding(false) to be a hadr./fragm.
  void SetParPossibility(G4bool choice);    // Set possibility(true)/forbiding(false) to be a parent
  void SetKMin(G4double kmin);              // Set k-minimal for the candidate
  void SetDenseProbability(G4double prep);    // Set dense-probability for the candidate 
  void SetPreProbability(G4double prep);    // Set pre-probability for the candidate 
  void SetRelProbability(G4double relP);    // Set the relative probility of hadronization
  void SetIntegProbability(G4double intP);  // Set integrated probability for randomization
  void SetSecondRelProb(G4double relP);     // Set 2nd relative probility of hadronization
  void SetSecondIntProb(G4double intP);     // Set 2nd integr. probability for randomization
  void SetMass(G4double newMass);           // Set 4Momentum at rest with the new Mass
										   
private:  								   
  G4bool   possible;                        // permission/forbiding preFlag to be a hadron/fragment
  G4bool   parPossible;                     // permission/forbiding preFlag to be a parent
  G4double kMin;                            // mu^2/2M (Q-case), ~BindingEnergy (Clust.-case)
  G4double denseProbability;                // a#of clusters of the type in dense region
  G4double preProbability;                  // a#of clusters of the type or Q-suppression
  G4double relativeProbability;             // relative probability of hadronization
  G4double integralProbability;             // integrated probability of randomization
  G4double secondRelProbability;            // spare relative probability of hadronization
  G4double secondIntProbability;            // spare integrated probability of randomization
  G4QParentClusterVector thePClusters;      // vector of parent clusters for candid.-fragment
};

inline G4int G4QCandidate::operator==(const G4QCandidate &right) const  {return this==&right;}
inline G4int G4QCandidate::operator!=(const G4QCandidate &right) const  {return this!=&right;}

inline G4QParentCluster* G4QCandidate::TakeParClust(G4int nPC)  {return thePClusters[nPC];}
inline G4int    G4QCandidate::GetPClustEntries()          const {return thePClusters.entries();}
inline G4bool   G4QCandidate::GetPossibility()            const {return possible;}
inline G4bool   G4QCandidate::GetParPossibility()         const {return parPossible;}
inline G4double G4QCandidate::GetKMin()                   const {return kMin;}
inline G4double G4QCandidate::GetDenseProbability()       const {return denseProbability;}
inline G4double G4QCandidate::GetPreProbability()         const {return preProbability;}
inline G4double G4QCandidate::GetRelProbability()         const {return relativeProbability;}
inline G4double G4QCandidate::GetIntegProbability()       const {return integralProbability;}
inline G4double G4QCandidate::GetSecondRelProb()          const {return secondRelProbability;}
inline G4double G4QCandidate::GetSecondIntProb()          const {return secondIntProbability;}

inline void G4QCandidate::ClearParClustVector()             {thePClusters.clearAndDestroy();}
inline void G4QCandidate::FillPClustVec(G4QParentCluster* pCl)    {thePClusters.insert(pCl);}
inline void G4QCandidate::SetPossibility(G4bool choice)           {possible=choice;}
inline void G4QCandidate::SetParPossibility(G4bool choice)        {parPossible=choice;}
inline void G4QCandidate::SetKMin(G4double kmin)                  {kMin=kmin;}
inline void G4QCandidate::SetDenseProbability(G4double prep)      {denseProbability=prep;}
inline void G4QCandidate::SetPreProbability(G4double prep)        {preProbability=prep;}
inline void G4QCandidate::SetRelProbability(G4double relP)        {relativeProbability=relP;}
inline void G4QCandidate::SetIntegProbability(G4double intP)      {integralProbability=intP;}
inline void G4QCandidate::SetSecondRelProb(G4double relP)         {secondRelProbability=relP;}
inline void G4QCandidate::SetSecondIntProb(G4double intP)         {secondIntProbability=intP;}
inline void G4QCandidate::SetMass(G4double newM)  {Set4Momentum(G4LorentzVector(0,0,0,newM));}

#endif
