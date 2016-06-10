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
// $Id: G4NuclWatcher.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 20100202  M. Kelsey -- Move most code into .cc file
// 20100405  M. Kelsey -- Pass const-ref std::vector<>
// 20101010  M. Kelsey -- Migrate to integer A and Z

#ifndef G4NUCL_WATCHER_HH
#define G4NUCL_WATCHER_HH

#include "G4Types.hh"

#include <algorithm>
#include <vector>
#include <cmath>

class G4NuclWatcher {
public:
  G4NuclWatcher(G4int z, 
		const std::vector<G4double>& expa, 
		const std::vector<G4double>& expcs, 
		const std::vector<G4double>& experr, 
		G4bool check, 
		G4bool nucl);

  ~G4NuclWatcher() {}

  void watch(G4int a, G4int z);
  void setInuclCs(G4double csec, G4int nev);

  G4double getChsq() const { return izotop_chsq; }
  G4bool to_check() const { return checkable; }
  G4bool look_forNuclei() const { return nucleable; }
  G4double getLhood() const { return aver_lhood; }
  G4double getNmatched() const { return aver_matched; }

  std::pair<G4double, G4double> getExpCs() const;
  std::pair<G4double, G4double> getInuclCs() const;

  std::pair<G4double, G4double> getAverageRatio() const { 
    return std::pair<G4double, G4double>(average_ratio, aver_rat_err); 
  }

  void print();

private: 
  G4int nuclz;
  G4double izotop_chsq;
  G4double average_ratio;
  G4double aver_rat_err;
  G4double aver_lhood;
  G4double aver_matched;
  std::vector<G4double> exper_as;
  std::vector<G4double> exper_cs;
  std::vector<G4double> exper_err;
  std::vector<G4double> simulated_as;
  std::vector<G4double> simulated_cs;
  std::vector<G4double> simulated_errors;
  std::vector<G4double> simulated_prob;
  G4bool checkable;
  G4bool nucleable;
};

#endif // G4NUCL_WATCHER_HH 

