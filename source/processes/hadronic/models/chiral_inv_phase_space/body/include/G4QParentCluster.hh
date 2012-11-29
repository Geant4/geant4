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
//      ---------------- G4QParentCluster ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for a Parent nuclear cluster in the CHIPS Model
// -------------------------------------------------------------------
// Short description: The parent cluster is the cluster, which can be
// used for the nuclear fragment production. Different clusters csn be
// used as the parent cluser for the particular G4QCandidate (nuclear
// fragment), e.g. t and He3 for the t-fragment production. So the
// G4QParentClusterVector is needed.
// -------------------------------------------------------------------

#ifndef G4QParentCluster_h
#define G4QParentCluster_h 1

#include <iostream>
#include "globals.hh"
#include "G4QContent.hh"

class G4QParentCluster
{
public:
  // Constructors
  G4QParentCluster(G4int PDGCode = 0);               // Construction by PDGCode
  G4QParentCluster(G4int PDGCode, G4double prob);    // Construction by PDGCode & Probab
  G4QParentCluster(const G4QParentCluster& rhs);     // Copy Constructor by value
  G4QParentCluster(G4QParentCluster* rhs);           // Copy Constructor by pointer

  ~G4QParentCluster();                               // Public Destructor

  // Operators
  const G4QParentCluster& operator=(const G4QParentCluster& rhs);
  G4bool                   operator==(const G4QParentCluster& rhs) const;
  G4bool                   operator!=(const G4QParentCluster& rhs) const;

  // Selectors
  G4int      GetPDGCode()     const; // Get PDG code of the Parent Cluster
  G4double   GetProbability() const; // Get a probability of hadronization on it
  G4int      GetNQPart2()     const; // Get n-2 for the fragment
  G4QContent GetTransQC()     const; // Get QuarkCont of a Pseudo Exchange Meson
  G4double   GetLow()         const; // Get a low limit for randomization
  G4double   GetHigh()        const; // Get a high limit for randomization
  G4double   GetEBMass()      const; // Get a Nuclear Bounded mass of the parent cluster
  G4double   GetEBind()       const; // Get Environment BindingEnergy for the parentCluster
  G4double   GetNBMass()      const; // Get Environmental BoundedMass of the parent Cluster
  G4double   GetNBind()       const; // Get TotalNucleusBindingEnergy for the parentCluster

  // Modifiers
  void  SetPDGCode(G4int newPDGCode);    // Set PDG code of the Parent Cluster
  void  SetProbability(G4double probab); // Set probab. to hadronize on this cluster
  void  SetNQPart2(G4int nm2);           // Get n-2 for the fragment
  void  SetTransQC(G4QContent newTrans); // Set QuarkCont of a Pseudo Exchange Meson
  void  SetLow(G4double loLim);          // Set a low limit for hadronization
  void  SetHigh(G4double hiLim);         // Set a high limit for hadronization
  void  SetEBMass(G4double bMass);       // Set a bounded mass of the parent cluster in E
  void  SetEBind(G4double bEn);          // Set binding energy of the parent cluster in E
  void  SetNBMass(G4double bMass);       // Set a bounded mass of the parent cluster in N
  void  SetNBind(G4double bEn);          // Set binding energy of the parent cluster in N

  // General

private:
  // Encapsulated functions

private:
  // the Body
  G4int               thePDGCode;
  G4double            theProbability;
  // Secondary
  G4int               nQPart2;
  G4QContent          transQC;           // Quark Content of the pseudo exchange meson
  G4double            lowLimit;
  G4double            highLimit;
  G4double            theEnvBoundedMass;
  G4double            theEnvBindingEnergy;
  G4double            theNucBoundedMass;
  G4double            theNucBindingEnergy;
};

// Not member operators
std::ostream&   operator<<(std::ostream& lhs, G4QParentCluster& rhs);
std::ostream&   operator<<(std::ostream& lhs, const G4QParentCluster& rhs);

inline G4bool G4QParentCluster::operator==(const G4QParentCluster& rhs) const
                                                                       {return this==&rhs;}
inline G4bool G4QParentCluster::operator!=(const G4QParentCluster& rhs) const
                                                                       {return this!=&rhs;}
 
inline G4int      G4QParentCluster::GetPDGCode()     const {return thePDGCode;}
inline G4double   G4QParentCluster::GetProbability() const {return theProbability;}
inline G4int      G4QParentCluster::GetNQPart2()     const {return nQPart2;}
inline G4QContent G4QParentCluster::GetTransQC()     const {return transQC;}
inline G4double   G4QParentCluster::GetHigh()        const {return highLimit;}
inline G4double   G4QParentCluster::GetLow()         const {return lowLimit;}
inline G4double   G4QParentCluster::GetEBMass()      const {return theEnvBoundedMass;}
inline G4double   G4QParentCluster::GetEBind()       const {return theEnvBindingEnergy;}
inline G4double   G4QParentCluster::GetNBMass()      const {return theNucBoundedMass;}
inline G4double   G4QParentCluster::GetNBind()       const {return theNucBindingEnergy;}

inline void  G4QParentCluster::SetPDGCode(G4int newPDGCode)    {thePDGCode   = newPDGCode;}
inline void  G4QParentCluster::SetProbability(G4double prob)   {theProbability = prob;}
inline void  G4QParentCluster::SetNQPart2(G4int nm2)           {nQPart2        = nm2;}
inline void  G4QParentCluster::SetTransQC(G4QContent newTrans) {transQC        = newTrans;}
inline void  G4QParentCluster::SetHigh(G4double hiLim)         {highLimit      = hiLim;}
inline void  G4QParentCluster::SetLow(G4double loLim)          {lowLimit       = loLim;}
inline void  G4QParentCluster::SetEBMass(G4double bMass)       {theEnvBoundedMass = bMass;}
inline void  G4QParentCluster::SetNBMass(G4double bMass)       {theNucBoundedMass = bMass;}
inline void  G4QParentCluster::SetEBind(G4double bEn)          {theEnvBindingEnergy= bEn;}
inline void  G4QParentCluster::SetNBind(G4double bEn)          {theNucBindingEnergy= bEn;}

#endif



