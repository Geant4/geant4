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
// $Id: G4MaterialCutsCouple.hh 70369 2013-05-29 14:59:24Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//
// Class Description
//  This class is 
//
// ------------------------------------------------------------
//   First Implementation          17 Sep. 2002  H.Kurahige
// ------------------------------------------------------------

#ifndef G4MaterialCutsCouple_h 
#define G4MaterialCutsCouple_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4Material;
#include "G4ProductionCuts.hh"

class G4MaterialCutsCouple  
{
  public: // with description
  //  constructor 
  G4MaterialCutsCouple();
  G4MaterialCutsCouple(const G4Material*, G4ProductionCuts* cut=0);

  //  copy constructor 
  G4MaterialCutsCouple(const G4MaterialCutsCouple &right);

  G4MaterialCutsCouple & operator=(const G4MaterialCutsCouple &right);

  public: 
  //  destructor 
  virtual ~G4MaterialCutsCouple();

  // equal opperators
  G4int operator==(const G4MaterialCutsCouple &right) const;
  G4int operator!=(const G4MaterialCutsCouple &right) const;

  public: // with description
  void              SetMaterial(const G4Material*);
  // Set pointer to material

  const G4Material* GetMaterial() const;
  // Get pointer to material

  void              SetProductionCuts(G4ProductionCuts*);
  // Set pointer to production cuts

  G4ProductionCuts* GetProductionCuts() const;
  // Get pointer to production cuts

  G4bool           IsRecalcNeeded() const;
  // return true if cut and/or material has been modified 
  // after last calculation of PhysicsTable          

  void             PhysicsTableUpdated();
  // inform end of calculation of PhysicsTable  

  void             SetIndex(G4int idx);
  G4int            GetIndex() const;
  // Set/Get the index number in G4ProductionCutsTable

  void             SetUseFlag(G4bool flg=true);
  G4bool           IsUsed() const;
 
  private:
  G4bool                   isMaterialModified;
  const G4Material*        fMaterial;
  G4ProductionCuts*        fCuts;
  G4int                    indexNumber;
  G4bool                   isUsedInGeometry;
};

#include "G4ProductionCuts.hh"
inline 
 void G4MaterialCutsCouple::SetIndex(G4int idx)
{ indexNumber = idx; }

inline
 G4int G4MaterialCutsCouple::GetIndex() const
{ return indexNumber; }

inline 
 void G4MaterialCutsCouple::SetUseFlag(G4bool flg)
{ isUsedInGeometry = flg; }

inline
 G4bool G4MaterialCutsCouple::IsUsed() const
{ return isUsedInGeometry; }

inline
 void G4MaterialCutsCouple::SetProductionCuts(G4ProductionCuts* aCut)
{ fCuts = aCut; }

inline
 G4ProductionCuts* G4MaterialCutsCouple::GetProductionCuts() const
{ return fCuts; }

inline
 G4int G4MaterialCutsCouple::operator==(const G4MaterialCutsCouple &right) const
{
  return (this == &right);
}

inline
 G4int G4MaterialCutsCouple::operator!=(const G4MaterialCutsCouple &right) const
{
  return (this !=  &right);
}

inline
  void  G4MaterialCutsCouple::SetMaterial(const G4Material* material)
{
  fMaterial = material;
  isMaterialModified = true;
}

inline
  const G4Material* G4MaterialCutsCouple::GetMaterial() const
{
  return fMaterial;
}

inline
G4bool  G4MaterialCutsCouple::IsRecalcNeeded() const
{
  G4bool isCutModified = false;
  if (fCuts !=0 ) isCutModified = fCuts->IsModified();
  return (isMaterialModified || isCutModified);
}

inline
void   G4MaterialCutsCouple::PhysicsTableUpdated()
{
  if (fCuts !=0 ) fCuts->PhysicsTableUpdated();
  isMaterialModified = false;
}


#endif






