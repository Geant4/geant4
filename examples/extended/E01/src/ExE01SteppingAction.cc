// $Id: ExE01SteppingAction.cc,v 1.2 1999/04/17 04:06:49 kurasige Exp $

#include "ExE01SteppingAction.hh"
#include "ExE01RunAction.hh"

void ExE01SteppingAction::UserSteppingAction(const G4Step*){

  HepRef(Histo1D) h1 = ExE01RunAction::get_1d();
  HepRef(Histo2D) h2 = ExE01RunAction::Get2d();

  HepRandom::setTheEngine(&theJamesEngine);

  double james = RandGauss::shoot(0.3,0.1);
  h1->fill(james,0.01);

  HepRandom::setTheEngine(&theDRand48Engine);

  double d48 = RandGauss::shoot(0.7,0.1);
  h2->fill(james,d48,0.01);   
}


