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
// $Id: G4ProductionCuts.hh 70369 2013-05-29 14:59:24Z gcosmo $
// GEANT4 tag $Name: geant4-09-04-ref-00 $
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
//   Add cuts for proton           28 Jul. 2009  H.Kurashige
// ------------------------------------------------------------

#ifndef G4ProductionCuts_h 
#define G4ProductionCuts_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>
#include "G4ParticleDefinition.hh"

enum G4ProductionCutsIndex
{
  idxG4GammaCut =0,
  idxG4ElectronCut,
  idxG4PositronCut,
 
  idxG4ProtonCut,  // for proton

  NumberOfG4CutIndex
};

class G4ProductionCuts  
{
  public: // with description
  //  constructor 
  G4ProductionCuts();

  //  copy constructor 
  G4ProductionCuts(const G4ProductionCuts &right);

  G4ProductionCuts & operator=(const G4ProductionCuts &right);

  public: 
  //  destructor 
  virtual ~G4ProductionCuts();

  // equal opperators
  G4int operator==(const G4ProductionCuts &right) const;
  G4int operator!=(const G4ProductionCuts &right) const;

  public: // with description
  // Set Cuts methods
  void              SetProductionCut(G4double cut, G4int index = -1);
  void              SetProductionCut(G4double cut, G4ParticleDefinition* ptcl);
  void              SetProductionCut(G4double cut, const G4String& pName);
  // Set the productionCut in range with an index to particle type
  // if index is omitted, the value is applied to all particles

  G4double          GetProductionCut(G4int index) const;
  // Get the productionCut in range with an index to particle type

  G4double          GetProductionCut(const G4String& name) const;
  // Get the productionCut in range with a name of particle type
  
  void              SetProductionCuts(std::vector<G4double>&);
  // Set the vector of production cuts in range for all particles

  const std::vector<G4double>&   GetProductionCuts() const;
  // Get the vector of production cuts in range for all particles

  G4bool           IsModified() const;
  // return true if any cut value has been modified 
  // after last calculation of PhysicsTable          

  void             PhysicsTableUpdated();
  // inform end of calculation of PhysicsTable  to ProductionCut 
 
  public:
  static G4int GetIndex(const G4String& name);
  static G4int GetIndex(const G4ParticleDefinition* ptcl);

  protected:
  std::vector<G4double>         fRangeCuts;
  G4bool                          isModified;

  private:
  static G4ThreadLocal G4ParticleDefinition* gammaDef;
  static G4ThreadLocal G4ParticleDefinition* electDef;
  static G4ThreadLocal G4ParticleDefinition* positDef;

  static G4ThreadLocal G4ParticleDefinition* protonDef; // for proton

};


inline
void  G4ProductionCuts::SetProductionCut(G4double cut, G4int index)
{
  if (index<0) {
    for(G4int i = 0; i < NumberOfG4CutIndex; i++) {
      fRangeCuts[i] = cut;
    }
    isModified = true;

  } else if (index < NumberOfG4CutIndex) {
    fRangeCuts[index] = cut;
    isModified = true;
  }     
}

inline 
void  G4ProductionCuts::SetProductionCut(G4double cut, G4ParticleDefinition* ptcl)
{
  G4int idx = -1;
  if(ptcl) idx = GetIndex(ptcl);
  if(idx>=0) SetProductionCut(cut,idx);
}

inline
void  G4ProductionCuts::SetProductionCut(G4double cut, const G4String& pName)
{
  G4int idx = GetIndex(pName);
  if(idx>=0) SetProductionCut(cut,idx);
}

inline
G4double  G4ProductionCuts::GetProductionCut(G4int index) const
{
  G4double cut=-1.0;
  if ( (index>=0) && (index<NumberOfG4CutIndex) ) {
    cut = fRangeCuts[index]; 
  }
  return cut;
}

inline
G4double  G4ProductionCuts::GetProductionCut(const G4String& name) const
{
  return GetProductionCut(GetIndex(name)); 
}


inline
const std::vector<G4double>&   G4ProductionCuts::GetProductionCuts() const
{
  return fRangeCuts;
}

inline
G4bool  G4ProductionCuts::IsModified() const
{
  return isModified;
}

inline
void   G4ProductionCuts::PhysicsTableUpdated()
{
  isModified = false;
}

#endif











