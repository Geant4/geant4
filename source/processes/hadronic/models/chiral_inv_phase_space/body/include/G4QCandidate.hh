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
// $Id$
//
//      ---------------- G4QCandidate ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for Quasmon initiated Candidates used by the CHIPS Model
// ----------------------------------------------------------------------
// Short description: A candidate for hadronization. The candidates
// (hadrons or nuclear fragments) are competative, each quark of a
// Quasmon select which candidate to use for hadronization
// ------------------------------------------------------------------

#ifndef G4QCandidate_h
#define G4QCandidate_h 1

#include "G4QHadron.hh"
#include "G4QParentClusterVector.hh"
#include <algorithm>

class G4QCandidate : public G4QHadron
{
public:
  G4QCandidate();                                            // Default Constructor
  G4QCandidate(G4int PDGcode);                               // Constructor by PDG Code
  G4QCandidate(const G4QCandidate& right);                   // Copy Constructor by value
  G4QCandidate(G4QCandidate* right);                         // Copy Constructor by pointer
  ~G4QCandidate();                                           // Public Destructor
  // Overloaded Operators
  const G4QCandidate& operator=(const G4QCandidate& right);
  G4bool operator==(const G4QCandidate &right) const;
  G4bool operator!=(const G4QCandidate &right) const;
  // Specific Selectors
  G4QParentCluster* TakeParClust(G4int nPC);// Get pointer to ParentClust from ParClastVect
  G4int    GetPClustEntries()    const;     // Get a#of Parent Clusters in ParClastVector
  G4bool   GetPossibility()      const;     // Get possibility(true)/forbid(false) hadr/fr
  G4bool   GetParPossibility()   const;     // Get possibility(true)/forbidi(false) parent
  G4double GetKMin()             const;     // Get k-minimal for the candidate
  G4double GetEBMass()           const;     // Get bound mass in respect to Environment 
  G4double GetNBMass()           const;     // Get bound mass in respect to TotalNucleus 
  G4double GetDenseProbability() const;     // Get dense-probability for the candidate 
  G4double GetPreProbability()   const;     // Get pre-probability for the candidate 
  G4double GetRelProbability()   const;     // Get the relative probility of hadronization
  G4double GetIntegProbability() const;     // Get integrated probability for randomization
  G4double GetSecondRelProb()    const;     // Get 2nd relative probility of hadronization
  G4double GetSecondIntProb()    const;     // Get 2nd integ. probability for randomization
  // Specific Modifiers
  void ClearParClustVector();               // Clear theParentClasterVector of theCandidate
  void FillPClustVec(G4QParentCluster* pCl);// Set pointer to ParentClust in ParClastVector
  void SetPossibility(G4bool choice);       // Set possibility(true)/forbid(false) hadr/fr
  void SetParPossibility(G4bool choice);    // Set possibility(true)/forbid(false) parent
  void SetKMin(G4double kmin);              // Set k-minimal for the candidate
  void SetDenseProbability(G4double prep);  // Set dense-probability for the candidate 
  void SetPreProbability(G4double prep);    // Set pre-probability for the candidate 
  void SetRelProbability(G4double relP);    // Set the relative probility of hadronization
  void SetIntegProbability(G4double intP);  // Set integrated probability for randomization
  void SetSecondRelProb(G4double relP);     // Set 2nd relative probility of hadronization
  void SetSecondIntProb(G4double intP);     // Set 2nd integrProbability for randomization
  void SetEBMass(G4double newMass);         // Set mass bounded to Environment
  void SetNBMass(G4double newMass);         // Set mass bounded to Total Nucleus

// Body             
private:
  G4bool   possible;                // permission/forbiding preFlag to be a hadron/fragment
  G4bool   parPossible;             // permission/forbiding preFlag to be a parent
  G4double kMin;                    // mu^2/2M (Q-case), ~BindingEnergy (Clust.-case)
  G4double denseProbability;        // a#of clusters of the type in dense region
  G4double preProbability;          // a#of clusters of the type or Q-suppression
  G4double relativeProbability;     // relative probability of hadronization
  G4double integralProbability;     // integrated probability of randomization
  G4double secondRelProbability;    // spare relative probability of hadronization
  G4double secondIntProbability;    // spare integrated probability of randomization
  G4QParentClusterVector thePClusters; // vector of parent clusters for candid.-fragment
  G4double EBMass;                     // Bound Mass of the cluster in the Environment
  G4double NBMass;                     // Bound Mass of the cluster in the Total Nucleus
};

inline G4bool G4QCandidate::operator==(const G4QCandidate &rhs) const  {return this==&rhs;}
inline G4bool G4QCandidate::operator!=(const G4QCandidate &rhs) const  {return this!=&rhs;}

inline G4QParentCluster* G4QCandidate::TakeParClust(G4int nPC){return thePClusters[nPC];}
inline G4int    G4QCandidate::GetPClustEntries()        const {return thePClusters.size();}
inline G4bool   G4QCandidate::GetPossibility()          const {return possible;}
inline G4bool   G4QCandidate::GetParPossibility()       const {return parPossible;}
inline G4double G4QCandidate::GetKMin()                 const {return kMin;}
inline G4double G4QCandidate::GetEBMass()               const {return EBMass;}
inline G4double G4QCandidate::GetNBMass()               const {return NBMass;}
inline G4double G4QCandidate::GetDenseProbability()     const {return denseProbability;}
inline G4double G4QCandidate::GetPreProbability()       const {return preProbability;}
inline G4double G4QCandidate::GetRelProbability()       const {return relativeProbability;}
inline G4double G4QCandidate::GetIntegProbability()     const {return integralProbability;}
inline G4double G4QCandidate::GetSecondRelProb()       const {return secondRelProbability;}
inline G4double G4QCandidate::GetSecondIntProb()       const {return secondIntProbability;}

inline void G4QCandidate::ClearParClustVector()                 
{
  std::for_each(thePClusters.begin(), thePClusters.end(), DeleteQParentCluster());
  thePClusters.clear();
}

inline void G4QCandidate::FillPClustVec(G4QParentCluster* pCl)
{
  thePClusters.push_back(pCl);                              // Fill new instance of PCl
}
inline void G4QCandidate::SetPossibility(G4bool choice)      {possible=choice;}
inline void G4QCandidate::SetParPossibility(G4bool choice)   {parPossible=choice;}
inline void G4QCandidate::SetKMin(G4double kmin)             {kMin=kmin;}
inline void G4QCandidate::SetDenseProbability(G4double prep) {denseProbability=prep;}
inline void G4QCandidate::SetPreProbability(G4double prep)   {preProbability=prep;}
inline void G4QCandidate::SetRelProbability(G4double relP)   {relativeProbability=relP;}
inline void G4QCandidate::SetIntegProbability(G4double intP) {integralProbability=intP;}
inline void G4QCandidate::SetSecondRelProb(G4double relP)    {secondRelProbability=relP;}
inline void G4QCandidate::SetSecondIntProb(G4double intP)    {secondIntProbability=intP;}
inline void G4QCandidate::SetEBMass(G4double newM)           {EBMass=newM;}
inline void G4QCandidate::SetNBMass(G4double newM)           {NBMass=newM;}

#endif
