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
// FermiBreakUp de-excitation model
// by V. Ivanchenko (July 2016)
//

#ifndef G4FermiDecayProbability_h
#define G4FermiDecayProbability_h 1

#include "globals.hh"

class G4FermiFragment;

class G4FermiDecayProbability 
{
public:

  explicit G4FermiDecayProbability();
  ~G4FermiDecayProbability() = default;

  G4double ComputeProbability(G4int Z, G4int A, G4int spin, G4double TotalE, 
			      const G4FermiFragment* f1, 
			      const G4FermiFragment* f2) const;

  G4FermiDecayProbability(const G4FermiDecayProbability &right) = delete;  
  const G4FermiDecayProbability & operator = 
  (const G4FermiDecayProbability &right) = delete;
 
};

#endif

