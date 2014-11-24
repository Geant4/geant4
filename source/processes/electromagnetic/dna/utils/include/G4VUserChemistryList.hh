// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

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

class G4Molecule;
class G4DNAMolecularReactionTable;
class G4VITStepModel;
class G4MoleculeDefinition;

class G4VPhysicsContructor;

class G4VUserChemistryList
{
public:
  G4VUserChemistryList();
  virtual ~G4VUserChemistryList();

  // If your user class also inherits from G4VPhysicsConstructor,
  // please put this flag to true
  virtual bool IsPhysicsConstructor()
  {
    if((G4VPhysicsContructor*)(this)) return true;
    return false;
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
  void RegisterTimeStepModel(G4VITStepModel* timeStepModel,
                             double startingTime = 0);
  void BuildPhysicsTable(G4MoleculeDefinition*);

  int verboseLevel;
};

#endif /* G4VUSERCHEMISTRYLIST_HH_ */
