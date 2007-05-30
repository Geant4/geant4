//******************************************************************************
// PrimaryGeneratorAction.hh
//
// This class is a class derived from G4VUserPrimaryGeneratorAction for 
// constructing the process used to generate incident particles.
//
// 1.00 JMV, LLNL, JAN-2007:  First version.
//******************************************************************************
// 
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "SponFissIsotope.hh"
#include "MultipleSource.hh"
#include "G4ThreeVector.hh"
#include "G4DataVector.hh"
#include "Randomize.hh"
#include "G4LLNLFission.hh"
#include "vector"

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    MultipleSource* fissionSource;    
    G4double time;
    G4double totalIntensity;
    G4int nisotopes;
};

#endif
