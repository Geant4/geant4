//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4MaterialCutsCouple.hh,v 1.2 2002-12-16 11:15:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
class G4ProductionCuts;

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
 
  private:
  G4bool                   isMaterialModified;
  const G4Material*        fMaterial;
  G4ProductionCuts*        fCuts;
};

#include "G4ProductionCuts.hh"

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
  return isMaterialModified || isCutModified;
}

inline
  void   G4MaterialCutsCouple::PhysicsTableUpdated()
{
  if (fCuts !=0 ) fCuts->PhysicsTableUpdated();
  isMaterialModified = true;
}

#endif






