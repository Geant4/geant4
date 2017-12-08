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
// $Id: G4MuonicAtom.hh 98732 2016-08-09 10:50:57Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      History: first implementation, 
//      July 2016, K. Lynch
//      June 2017, K.L. Genser added baseion, lifetimes and access functions

#ifndef G4MuonicAtom_h 
#define G4MuonicAtom_h 1

#include "G4Ions.hh"

class G4MuonicAtom : public G4Ions
{
 protected:
   G4MuonicAtom(){};

 public:
  G4MuonicAtom(
               const G4String&     aName,        G4double            mass,
               G4double            width,        G4double            charge,   
               G4int               iSpin,        G4int               iParity,    
               G4int               iConjugation, G4int               iIsospin,   
               G4int               iIsospin3,    G4int               gParity,
               const G4String&     pType,        G4int               lepton,      
               G4int               baryon,       G4int               encoding,
               G4bool              stable,       G4double            lifeTime,
               G4DecayTable        *decaytable,  G4bool              shortlived,
               const G4String&     subType,      G4Ions const*       baseion,
               G4int               anti_encoding  =0,
               G4double            excitation = 0.0, 
               G4int               isomer = 0,
               G4double            DIOLifeTime = -1.0,
               G4double            NCLifeTime = -1.0
               );

  virtual                      ~G4MuonicAtom();
  G4MuonicAtom*                MuonicAtomDefinition();
  G4MuonicAtom*                MuonicAtom();
  G4Ions const*                GetBaseIon() const;
  G4double                     GetDIOLifeTime() const;
  void                         SetDIOLifeTime(G4double lt);
  G4double                     GetNCLifeTime() const;
  void                         SetNCLifeTime(G4double lt);

 private:
   G4Ions const*                baseIon;
   G4double                     fDIOLifeTime;
   G4double                     fNCLifeTime; 
};

inline
 G4MuonicAtom* G4MuonicAtom::MuonicAtomDefinition()
{
  return this;
}

inline
 G4MuonicAtom* G4MuonicAtom::MuonicAtom() 
{
  return this;
}

inline
 G4Ions const* G4MuonicAtom::GetBaseIon() const
{
  return baseIon;
}

inline
 G4double G4MuonicAtom::GetDIOLifeTime() const
{
  return fDIOLifeTime;
}

inline
 void G4MuonicAtom::SetDIOLifeTime(G4double lt)
{
  fDIOLifeTime = lt;
}

inline
 G4double G4MuonicAtom::GetNCLifeTime() const
{
  return fNCLifeTime;
}

inline
 void G4MuonicAtom::SetNCLifeTime(G4double lt)
{
  fNCLifeTime = lt;
}

#endif
