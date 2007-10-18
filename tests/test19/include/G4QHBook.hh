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

#ifndef G4QHBook_h
#define G4QHBook_h 1

// ----------------------------------------------------------------------
//      GEANT 4 class header file
//      ---------------- G4QHBook ----------------
//      made for process level tests by Mikhail Kossov - Feb 2005
//      class for booking and filling histograms and ntuples
// ----------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "G4VParticleChange.hh"
#include "G4DynamicParticle.hh"

class G4QHBook
{
public:
  G4QHBook();           // Default Constructor
  ~G4QHBook();          // Destructor

  // Specific Modifiers
  void  FillEvt(const G4VParticleChange* hadrons, const G4DynamicParticle* pdProj);

private: 
  G4int         nEvnt;    // Consecutive number of call to fill the histograms 
  std::ofstream histNevt; // 1D histogram, id=1, to store number of calls to fill
  std::ofstream tuplEvtA; // Ntuple, id=25, to be filled once per event (All particles)
  std::ofstream tuplEvtQ; // Ntuple, id=27, to be filled once per event (Quasmon particles)
  std::ofstream tuplIncl; // Ntuple, id=20, to be filled once per track
  std::ofstream tuple3pi; // Ntuple, id=22, to be filled once per event with 3 pions
};

#endif
