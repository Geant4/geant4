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
// $Id: G4ParticleWithCuts.cc,v 1.18 2002-12-16 11:15:45 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//  History: 
//   first implementation, based on object model of Hisaya Kurashige,
//                                                  21 Oct 1996
//   calculation of Range Table is based on implementeation for Muon 
//                                         by L.Urban, 10 May 1996
//   modify CalcEnergyCuts                 09 Nov. 1998, L.Urban
//   added  RestoreCuts  H.Kurashige 09 Mar. 2001
//   modify for material-V03-02-02 (STL migration)  H.Kurashige 19 Sep. 2001
//   introduced material dependent range cuts   08 Oct. 2001
// ------------------------------------------------------------
#include "globals.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"

#include "G4ProductionCutsTable.hh"
#include "G4VRangeToEnergyConverter.hh"

#include "G4ios.hh"

#include "g4std/strstream"

G4ParticleWithCuts::G4ParticleWithCuts(
		const G4String&  aName,  
                G4double         mass,     
                G4double         width,
                G4double         charge,   
                G4int            iSpin,
                G4int            iParity,
                G4int            iConjugation,
                G4int            iIsospin,   
                G4int            iIsospinZ, 
                G4int            gParity,
                const G4String&  pType,
                G4int            lepton,
                G4int            baryon,
                G4int            encoding,
                G4bool           stable,
                G4double         lifetime,
                G4DecayTable     *decaytable,
		G4bool           shortlived)
	: G4ParticleDefinition(aName, mass, width, charge, iSpin, iParity,
                               iConjugation, iIsospin, iIsospinZ, gParity,
                               pType, lepton, baryon, encoding, stable,
                               lifetime, decaytable, shortlived), 
          theCutInMaxInteractionLength(0),
          theKineticEnergyCuts(0)
{
  // -- set ApplyCutsFlag in default   ----------
  SetApplyCutsFlag(false);
  
  // -- Production Cuts Table --
  theCutsTable =   G4ProductionCutsTable::GetProductionCutsTable();
  
}

G4ParticleWithCuts::~G4ParticleWithCuts()
{ 
  if(theCutInMaxInteractionLength !=0) delete [] theCutInMaxInteractionLength;
  if(theKineticEnergyCuts !=0) delete [] theKineticEnergyCuts;
}

// **********************************************************************
// **************************** SetCuts *********************************
// **********************************************************************

G4int  G4ParticleWithCuts::GetParticleIndex() const
{ 
  G4String name = GetParticleName();
  G4int index;
  if       ( name == "gamma" )        { index =  0; }
  else  if ( name == "e-" )           { index =  1; }
  else  if ( name == "e+" )           { index =  2; }
  else  if ( name == "proton" )       { index =  3; }
  else  if ( name == "anti_proton" )  { index =  4; }
  else  if ( name == "neutron" )      { index =  5; }
  else  if ( name == "anti_neutron" ) { index =  6; }
  else                               { index = -1; }
  return index;
}

void G4ParticleWithCuts::SetCuts(G4double aCut)
{
  G4int index = GetParticleIndex(); 
  if (index<0) {
//  #ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << " G4ParticleWithCuts::SetCuts ";
      G4cout << " Production Cut is not defined for ";
      G4cout << "[" << GetParticleName() <<"]" << G4endl;
    }
//  #endif
    return;
  }

  // Get Production Cuts for WORLD     
  // const G4MaterialCutsCouple* aCouple = theCutsTable->GetMaterialCutsCouple(0);
  // G4ProductionCuts* defaultCuts = aCouple->GetProductionCuts();
  G4ProductionCuts* defaultCuts
    = G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts();

  // Set cut in range to ProductionCuts  
  defaultCuts->SetProductionCut(aCut, index);

  return;
}

void G4ParticleWithCuts::SetRangeCut(G4double aCut, const G4Material* aMaterial)
{
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << " G4ParticleWithCuts::SetRangeCuts is obsolete !!" << G4endl;
      G4cout << " Use Production Cuts for Regions " << G4endl;
    }
#endif
}

void G4ParticleWithCuts::SetRangeCutVector(G4std::vector<G4double>& cuts)
{
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << " G4ParticleWithCuts::SetRangeCutVector is obsolete !!" << G4endl;
      G4cout << " Use Production Cuts for Regions " << G4endl;
    }
#endif 
}    

void G4ParticleWithCuts::SetEnergyRange(G4double lowedge, G4double highedge)
{
  G4VRangeToEnergyConverter::SetEnergyRange(lowedge, highedge);    
}


G4double* G4ParticleWithCuts::GetLengthCuts() const
{
  G4int idx = GetParticleIndex(); 
  if (idx<0) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << " G4ParticleWithCuts::GetLengthCuts ";
      G4cout << " Production Cut is not defined for ";
      G4cout << "[" << GetParticleName() <<"]" << G4endl;
    }
#endif
    return 0;
  }

  return G4ProductionCutsTable::GetProductionCutsTable()
              ->GetRangeCutsDoubleVector(size_t(idx));
}


G4double* G4ParticleWithCuts::GetEnergyCuts() const
{
   G4int idx = GetParticleIndex(); 
  if (idx<0) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << " G4ParticleWithCuts::GetEnergyCuts ";
      G4cout << " Production Cut is not defined for ";
      G4cout << "[" << GetParticleName() <<"]" << G4endl;
    }
#endif
    return 0;
  }
  
  return  G4ProductionCutsTable::GetProductionCutsTable()
              ->GetEnergyCutsDoubleVector(size_t(idx));
}
  

void G4ParticleWithCuts::ResetCuts()
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
    G4cout << " G4ParticleWithCuts::ResetCuts is obsolete !!" << G4endl;
    G4cout << " Use Production Cuts for Regions " << G4endl;
  }
#endif 
}


void G4ParticleWithCuts::ReCalcCuts()
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
    G4cout << " G4ParticleWithCuts::ResetCalcCuts is obsolete !!" << G4endl;
    G4cout << " Use Production Cuts for Regions " << G4endl;
  }
#endif 
}

G4double G4ParticleWithCuts::GetEnergyThreshold(const G4Material* aMaterial) const
{
   G4int index = GetParticleIndex(); 
  if (index<0) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << " G4ParticleWithCuts::GetEnergyThreshold ";
      G4cout << " Production Cut is not defined for ";
      G4cout << "[" << GetParticleName() <<"]" << G4endl;
    }
#endif
    return -1;
  }
  
  const G4std::vector<G4double>* energyVector = theCutsTable->GetEnergyCutsVector(index);
  return (*energyVector)[aMaterial->GetIndex()]; 
}


G4double G4ParticleWithCuts::GetRangeThreshold(const G4Material* aMaterial) const
{
   G4int index = GetParticleIndex(); 
  if (index<0) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << " G4ParticleWithCuts::GetRangeThreshold ";
      G4cout << " Production Cut is not defined for ";
      G4cout << "[" << GetParticleName() <<"]" << G4endl;
    }
#endif
    return -1;
  }
  
  const G4std::vector<G4double>* rangeVector = theCutsTable->GetRangeCutsVector(index);
  return (*rangeVector)[aMaterial->GetIndex()]; 
}
