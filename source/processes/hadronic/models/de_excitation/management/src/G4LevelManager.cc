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
// $Id: G4LevelManager.cc 88375 2015-02-16 17:31:21Z vnivanch $
//
// -------------------------------------------------------------------
//
//      GEANT4 source file 
//
//      File name:     G4LevelManager
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 4 January 2012
//
//      Modifications: 
//  13.02.2015 Design change for gamma de-excitation 
//      
// -------------------------------------------------------------------

#include "G4LevelManager.hh"
#include "G4HadronicException.hh"

G4String G4LevelManager::fFloatingLevels[] = {
  "-", "+X", "+Y", "+Z", "+U", "+V", "+W", "+R", "+S", "+T", "+A", "+B", "+C"};

G4LevelManager::G4LevelManager(size_t ntrans,
                               const std::vector<G4double>& energies,
			       const std::vector<G4int>& spin,
			       const std::vector<const G4NucLevel*>& levels)
  : nTransitions(0)
{ 
  if(0 < ntrans) {
    nTransitions = ntrans - 1;
    fLevelEnergy.reserve(ntrans);
    fSpin.reserve(ntrans);
    fLevels.reserve(ntrans);
    for(size_t i=0; i<ntrans; ++i) {
      fLevelEnergy.push_back(energies[i]);
      fSpin.push_back(spin[i]);
      fLevels.push_back(levels[i]);
    }
    //G4cout << "New G4LevelManager N= " << nTransitions << " " 
    //<< fLevelEnergy.size() << " <" << this << ">" << G4endl;
  }
}

G4LevelManager::~G4LevelManager()
{
  for(size_t i=0; i<=nTransitions; ++i) { delete fLevels[i]; }
}

size_t  
G4LevelManager::NearestLevelIndex(G4double energy, size_t index) const
{
  //G4cout<< "index= " << index << " max= " << nTransitions << " exc= " << ener 
  //	 << " Emax= " << fLevelEnergy[nTransitions] << G4endl;
  size_t idx = std::min(index, nTransitions);
  static const G4double tolerance = 1.0f-6;
  if(0 == nTransitions || std::abs(energy - fLevelEnergy[idx]) <= tolerance) {
    return idx;
  }
  // ground state
  if(energy <= fLevelEnergy[1]*0.5)  
    { idx = 0; }
  // take top level
  else if((fLevelEnergy[nTransitions] + fLevelEnergy[nTransitions-1])*0.5 <= energy) 
    { idx = nTransitions; }

  // if shortcuts are not working, make binary search
  else {
    idx = std::lower_bound(fLevelEnergy.begin(), fLevelEnergy.end(), energy)
      - fLevelEnergy.begin() - 1;
    if(energy - fLevelEnergy[idx] > fLevelEnergy[idx+1] - energy) { ++idx; }
    //G4cout << "E= " << energy << " " << fLevelEnergy[idx-1] 
    //<< " " << fLevelEnergy[idx] << G4endl;
  }
  return idx;
}

const G4String& G4LevelManager::FloatingType(size_t i) const
{
#ifdef G4VERBOSE
  if(i > nTransitions) { PrintError(i, "Meta"); }
#endif
  return fFloatingLevels[fSpin[i]/100000]; 
  //  return fFloatingLevels[fFloating[i]]; 
}

#ifdef G4VERBOSE
void G4LevelManager::PrintError(size_t idx, const G4String& ss) const
{
  G4String sss = "G4LevelManager::"+ss+"()";
  G4ExceptionDescription ed;
  ed << "Index of a level " << idx << " > " 
     << nTransitions << " (number of levels)";
  G4Exception(sss,"had061",JustWarning,ed,"");
  throw G4HadronicException(__FILE__, __LINE__,"FATAL Hadronic Exception");
}  
#endif
