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
// $Id: G4ProductionCuts.hh,v 1.2 2002-12-16 11:15:43 gcosmo Exp $
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

#ifndef G4ProductionCuts_h 
#define G4ProductionCuts_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/vector"

enum G4ProductionCutsIndex
{
  idxG4GammaCut =0,
  idxG4ElectronCut,
  idxG4PositronCut,
  idxG4ProtonCut,
  idxG4AntiProtonCut,
  idxG4NeutronCut,
  idxG4AntiNeutronCut,

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
  // Set the productionCut in range with an index to particle type
  // if index is omitted, the value is applied to all particles

  G4double          GetProductionCut(G4int index) const;
  // Get the productionCut in range with an index to particle type

  G4double          GetProductionCut(const G4String& name) const;
  // Get the productionCut in range with a name of particle type
  
  void              SetProductionCuts(G4std::vector<G4double>&);
  // Set the vector of production cuts in range for all particles

  const G4std::vector<G4double>&   GetProductionCuts() const;
  // Get the vector of production cuts in range for all particles

  G4bool           IsModified() const;
  // return true if any cut value has been modified 
  // after last calculation of PhysicsTable          

  void             PhysicsTableUpdated();
  // inform end of calculation of PhysicsTable  to ProductionCut 
 
  protected:
  G4int            GetIndex(const G4String& name) const;

  protected:
  G4std::vector<G4double>         fRangeCuts;
  G4bool                          isModified;
};

#include "G4ProductionCuts.icc"

#endif











