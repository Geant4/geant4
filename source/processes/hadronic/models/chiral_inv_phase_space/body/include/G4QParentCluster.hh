// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QParentCluster.hh,v 1.2 2000-09-04 07:44:01 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4QParentCluster_h
#define G4QParentCluster_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QParentCluster ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for a Parent nuclear cluster in the CHIPS Model
// -------------------------------------------------------------

#include <iostream.h>
#include "globals.hh"
#include "G4QContent.hh"

class G4QParentCluster
{
public:
  // Constructors
  G4QParentCluster(G4int PDGCode = 0);               // Construction by PDGCode
  G4QParentCluster(G4int PDGCode, G4double prob);    // Construction by PDGCode & Probab
  G4QParentCluster(const G4QParentCluster& rhs);     // Copy Constructor

  ~G4QParentCluster();                               // Destructor

  // Operators
  const G4QParentCluster& operator=(const G4QParentCluster& rhs);
  G4int                   operator==(const G4QParentCluster& rhs) const;
  G4int                   operator!=(const G4QParentCluster& rhs) const;

  // Selectors
  G4int      GetPDGCode()      const;   // Get PDG code of the Parent Cluster
  G4double   GetProbability()  const;   // Get a probability of hadronization on it
  G4int      GetNQPart2()      const;   // Get n-2 for the fragment
  G4QContent GetTransQC()      const;   // Get QuarkCont of a Pseudo Exchange Meson
  G4double   GetLow()          const;   // Get a low limit for randomization
  G4double   GetHigh()         const;   // Get a high limit for randomization
  G4double   GetMass()         const;   // Get a bounded mass of the parent cluster
  G4double   GetBind()         const;   // Get binding energy of the parent cluster

  // Modifiers
  void  SetPDGCode(G4int newPDGCode);                // Set PDG code of the Parent Cluster
  void  SetProbability(G4double probab);             // Set probab. to hadronize on this cluster
  void  SetNQPart2(G4int nm2);                       // Get n-2 for the fragment
  void  SetTransQC(G4QContent newTrans);             // Set QuarkCont of a Pseudo Exchange Meson
  void  SetLow(G4double loLim);                      // Set a low limit for hadronization
  void  SetHigh(G4double hiLim);                     // Set a high limit for hadronization
  void  SetMass(G4double bMass);                     // Set a bounded mass of the parent cluster
  void  SetBind(G4double bEn);                       // Set binding energy of the parent cluster

  // General

private:
  // Encapsulated functions

private:
  // the Body
  G4int               thePDGCode;
  G4double            theProbability;
  // Secondary
  G4int               nQPart2;
  G4QContent          transQC;                         // Quark Content of pseudo exchange meson
  G4double            lowLimit;
  G4double            highLimit;
  G4double            theBoundedMass;
  G4double            theBindingEnergy;
};

// Not member operators
ostream&   operator<<(ostream& lhs, G4QParentCluster& rhs);
ostream&   operator<<(ostream& lhs, const G4QParentCluster& rhs);

inline G4int G4QParentCluster::operator==(const G4QParentCluster& rhs) const {return this==&rhs;}
inline G4int G4QParentCluster::operator!=(const G4QParentCluster& rhs) const {return this!=&rhs;}
 
inline G4int      G4QParentCluster::GetPDGCode()     const {return thePDGCode;}
inline G4double   G4QParentCluster::GetProbability() const {return theProbability;}
inline G4int      G4QParentCluster::GetNQPart2()     const {return nQPart2;}
inline G4QContent G4QParentCluster::GetTransQC()     const {return transQC;}
inline G4double   G4QParentCluster::GetHigh()        const {return highLimit;}
inline G4double   G4QParentCluster::GetLow()         const {return lowLimit;}
inline G4double   G4QParentCluster::GetMass()        const {return theBoundedMass;}
inline G4double   G4QParentCluster::GetBind()        const {return theBindingEnergy;}

inline void  G4QParentCluster::SetPDGCode(G4int newPDGCode)    {thePDGCode      = newPDGCode;}
inline void  G4QParentCluster::SetProbability(G4double prob)   {theProbability  = prob;}
inline void  G4QParentCluster::SetNQPart2(G4int nm2)           {nQPart2         = nm2;}
inline void  G4QParentCluster::SetTransQC(G4QContent newTrans) {transQC=newTrans;}
inline void  G4QParentCluster::SetHigh(G4double hiLim)         {highLimit       = hiLim;}
inline void  G4QParentCluster::SetLow(G4double loLim)          {lowLimit        = loLim;}
inline void  G4QParentCluster::SetMass(G4double bMass)         {theBoundedMass  = bMass;}
inline void  G4QParentCluster::SetBind(G4double bEn)           {theBindingEnergy= bEn;}

#endif



