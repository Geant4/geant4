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
// Created 25.03.2025 V.Ivanchenko 
// on base of codes of S.Incerti & M.Karamitros
//
// Double differential cross section data structure
//

#ifndef  G4DNASamplingTable_HH
#define  G4DNASamplingTable_HH 1

#include "globals.hh"
#include <vector>

class G4DNASamplingTable
{ 
public:
  explicit G4DNASamplingTable(std::size_t npoint);
  ~G4DNASamplingTable();

  void LoadData(const G4String& filename, G4double factorE, G4double scaleFactor,
		G4bool verbose);

  G4double GetValue(G4double ekinPrimary, G4double ekinSecondary, G4int shell) const;

  G4double SampleCumulative(G4double ekinPrimary, G4int shell) const;
    
  G4DNASamplingTable(const G4DNASamplingTable & copy) = delete;
  G4DNASamplingTable& operator=(const G4DNASamplingTable& right) = delete;

private:

  G4int GetIndex(const std::vector<G4double>&, G4double x) const;
  
  G4double VecInterpolation(const std::vector<G4double>* ener,
			    const std::vector<G4double>* val, G4double energy) const;
  
  G4double Interpolate(G4double e1, G4double e2, G4double e,
		       G4double xs1, G4double xs2) const;

  G4int fNpoints{0};
  std::vector<G4double> fPrimaryEnergy;
  std::vector<std::vector<G4double>* > fSecEnergy;
  std::vector<std::vector<G4double>* > fPDF[5];
  
};
#endif
