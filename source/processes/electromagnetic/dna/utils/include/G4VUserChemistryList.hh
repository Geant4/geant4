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
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4VUSERCHEMISTRYLIST_HH_
#define G4VUSERCHEMISTRYLIST_HH_

#include "G4Types.hh"
class G4Molecule;
class G4DNAMolecularReactionTable;
class G4VITStepModel;
class G4MoleculeDefinition;

enum TimeStepModel
{
  fSBS,
  fIRT,
  fIRT_syn
};

class G4VUserChemistryList
{
public:
  G4VUserChemistryList(G4bool flag = true);
  virtual ~G4VUserChemistryList();

  // If your user class also inherits from G4VPhysicsConstructor,
  // please put this flag to true
  bool IsPhysicsConstructor()
  {
    return fIsPhysicsConstructor;
  }
  
  void ThisIsAPhysicsConstructor(G4bool flag = true)
  {
    fIsPhysicsConstructor = flag;
  }

  ////////////////////////////////
  // to be called from PhysicsList

  virtual void ConstructMolecule()
  {
    ;
  } // PhysicsList::ConstructParticle
  virtual void ConstructProcess()
  {
    ;
  } // PhysicsList::ConstructProcess

  /////////////

  virtual void ConstructDissociationChannels()
  {
    ;
  }
  virtual void ConstructReactionTable(G4DNAMolecularReactionTable* reactionTable) = 0;
  virtual void ConstructTimeStepModel(G4DNAMolecularReactionTable* reactionTable) = 0;

  void BuildPhysicsTable();

 protected:
  void BuildPhysicsTable(G4MoleculeDefinition*);

  void RegisterTimeStepModel(G4VITStepModel* timeStepModel,
                             G4double startingTime = 0);

  G4int verboseLevel;
  G4bool fIsPhysicsConstructor;
};

#endif /* G4VUSERCHEMISTRYLIST_HH_ */
