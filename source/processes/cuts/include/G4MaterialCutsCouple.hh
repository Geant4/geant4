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
// G4MaterialCutsCouple
//
// Class description:
//
// The same material can be used in regions with different cut values,
// and physics processes must prepare several different cross-sections
// for that material.
// The G4ProductionCutsTable has G4MaterialCutsCouple objects, each of
// which consists of a material paired with a cut value.
// G4MaterialCutsCouples objects are numbered with an index which is the
// same as the index of a G4PhysicsVector for the corresponding
// G4MaterialCutsCouple in the G4PhysicsTable.
// The list of G4MaterialCutsCouple objects used in the geometry setup
// is updated before starting the event loop in each run.

// Author: H.Kurashige, 17 September 2002 - First implementation
// --------------------------------------------------------------------
#ifndef G4MaterialCutsCouple_hh 
#define G4MaterialCutsCouple_hh 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ProductionCuts.hh"

class G4Material;

class G4MaterialCutsCouple  
{
  public:

    G4MaterialCutsCouple();
    G4MaterialCutsCouple(const G4Material*, G4ProductionCuts* cut = nullptr);
      // Constructors

    virtual ~G4MaterialCutsCouple();
      // Destructor 

    G4MaterialCutsCouple(const G4MaterialCutsCouple& right);
    G4MaterialCutsCouple& operator=(const G4MaterialCutsCouple& right);
      // Copy constructor & assignment operator

    G4bool operator==(const G4MaterialCutsCouple& right) const;
    G4bool operator!=(const G4MaterialCutsCouple& right) const;
      // Equality operators

    void SetMaterial(const G4Material*);
      // Set pointer to material

    const G4Material* GetMaterial() const;
      // Get pointer to material

    void SetProductionCuts(G4ProductionCuts*);
      // Set pointer to production cuts

    G4ProductionCuts* GetProductionCuts() const;
      // Get pointer to production cuts

    G4bool IsRecalcNeeded() const;
      // Return true if cut and/or material has been modified 
      // after last calculation of PhysicsTable          

    void PhysicsTableUpdated();
      // Inform end of calculation of Physics Table  

    void SetIndex(G4int idx);
    G4int GetIndex() const;
      // Set/Get the index number in G4ProductionCutsTable

    void SetUseFlag(G4bool flg = true);
    G4bool IsUsed() const;
 
  private:

    G4bool            isMaterialModified = false;
    const G4Material* fMaterial = nullptr;
    G4ProductionCuts* fCuts = nullptr;
    G4int             indexNumber = -1;
    G4bool            isUsedInGeometry = false;
};

// ------------------
// Inline methods
// ------------------

inline 
void G4MaterialCutsCouple::SetIndex(G4int idx)
{
  indexNumber = idx;
}

inline
G4int G4MaterialCutsCouple::GetIndex() const
{
  return indexNumber;
}

inline 
void G4MaterialCutsCouple::SetUseFlag(G4bool flg)
{
  isUsedInGeometry = flg;
}

inline
G4bool G4MaterialCutsCouple::IsUsed() const
{
  return isUsedInGeometry;
}

inline
void G4MaterialCutsCouple::SetProductionCuts(G4ProductionCuts* aCut)
{
  fCuts = aCut;
}

inline
G4ProductionCuts* G4MaterialCutsCouple::GetProductionCuts() const
{
  return fCuts;
}

inline
G4bool G4MaterialCutsCouple::operator==(const G4MaterialCutsCouple& right) const
{
  return (this == &right);
}

inline
G4bool G4MaterialCutsCouple::operator!=(const G4MaterialCutsCouple& right) const
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
  if (fCuts != nullptr ) isCutModified = fCuts->IsModified();
  return (isMaterialModified || isCutModified);
}

inline
void G4MaterialCutsCouple::PhysicsTableUpdated()
{
  if (fCuts != nullptr ) fCuts->PhysicsTableUpdated();
  isMaterialModified = false;
}

#endif
