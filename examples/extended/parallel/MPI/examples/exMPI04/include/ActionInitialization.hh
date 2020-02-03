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
/// @file ActionInitialization.hh 
/// @brief Define action initialization

#ifndef ACTION_INITIALIZATION_H
#define ACTION_INITIALIZATION_H

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

<<<<<<< HEAD:source/processes/hadronic/models/parton_string/qgsm/include/G4GammaAnnCrossSection.hh
class G4GammaAnnCrossSection : public G4VAnnihilationCrossSection
{
public:
	G4GammaAnnCrossSection();
	G4bool InCharge(G4int aCode, G4int bCode);
	G4double GetXsec(G4double S);
	virtual ~G4GammaAnnCrossSection(){}

private:

	std::vector<G4ASCCrossSection*> theGammaNucXSections;
=======
class ActionInitialization : public G4VUserActionInitialization {
public:
  ActionInitialization(G4bool useNtuple, G4bool mergeNtuple);
  ~ActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;
private:
  G4bool fUseNtuple;  
  G4bool fMergeNtuple;  
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/parallel/MPI/examples/exMPI04/include/ActionInitialization.hh
};

#endif
