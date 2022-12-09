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
#include "G4NuclearLevelData.hh"
#include "G4ShellCorrection.hh"
#include "G4HadronicException.hh"
#include "G4Pow.hh"
#include <iomanip>
#include <CLHEP/Units/SystemOfUnits.h>

G4String G4LevelManager::fFloatingLevels[] = {
  "-", "+X", "+Y", "+Z", "+U", "+V", "+W", "+R", "+S", "+T", "+A", "+B", "+C"};

G4LevelManager::G4LevelManager(G4int Z, G4int A, std::size_t ntrans,
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
    for(std::size_t i=0; i<ntrans; ++i) {
      fLevelEnergy.push_back(energies[i]);
      fSpin.push_back(spin[i]);
      fLevels.push_back(levels[i]);
    }
    //G4cout << "New G4LevelManager N= " << nTransitions << " " 
    //<< fLevelEnergy.size() << " <" << this << ">" << G4endl;
  }
  auto ndata = G4NuclearLevelData::GetInstance();
  fLevelDensity = ndata->GetLevelDensity(Z, A, 0.0);

  // J. Nucl. Sci. Tech. 31(2): 151-162 (1994)
  fShellCorrection = ndata->GetShellCorrection()->GetShellCorrection(A,Z);
  if(A > 20) {
    G4int N = A - Z;
    G4int In = N - (N/2)*2; 
    G4int Iz = Z - (Z/2)*2;
    G4double a13 = 1.0/G4Pow::GetInstance()->Z13(A);
    if(In == 0 && Iz == 0) {
      fLevelDensity = 0.067946*A*(1.0 + 4.1277*a13);
    } else if(In == 0 && Iz == 1) {
      fLevelDensity = 0.053061*A*(1.0 + 7.1862*a13);
    } else if(In == 1 && Iz == 0) {
      fLevelDensity = 0.060920*A*(1.0 + 3.8767*a13);
    } else {
      fLevelDensity = 0.065291*A*(1.0 + 4.4505*a13);
    }
  }
}

G4LevelManager::~G4LevelManager()
{
  for(std::size_t i=0; i<=nTransitions; ++i) { delete fLevels[i]; }
}

std::size_t G4LevelManager::NearestLevelIndex(const G4double energy,
                                              const std::size_t index) const
{
  std::size_t idx = std::min(index, nTransitions);
  static const G4double tolerance = 10*CLHEP::eV;
  if(0 == nTransitions || std::abs(energy - fLevelEnergy[idx]) <= tolerance) {
    return idx;
  }
  idx = NearestLowEdgeLevelIndex(energy);
  if(idx < nTransitions && 
     (fLevelEnergy[idx] + fLevelEnergy[idx+1])*0.5 <= energy) { ++idx; }

  return idx;
}

const G4String& G4LevelManager::FloatingType(const std::size_t i) const
{
  return fFloatingLevels[fSpin[i]/100000]; 
}

void G4LevelManager::StreamInfo(std::ostream& out) const
{
  for(std::size_t i=0; i<=nTransitions; ++i) {
    G4long prec = out.precision(6);
    out << std::setw(6) << i << ". " 
	<< std::setw(8) << fLevelEnergy[i];
    if(fLevels[i]) {
	out << std::setw(8) << fLevels[i]->GetTimeGamma()
	    << std::setw(4) << fLevels[i]->NumberOfTransitions()
	    << std::setw(4) << SpinTwo(i)
	    << std::setw(4) << Parity(i)
	    << std::setw(4) << FloatingLevel(i);
    }
    out << "\n";
    out.precision(prec);
    if(fLevels[i]) { fLevels[i]->StreamInfo(out); }
  }
}
