#include "MyRunAction.hh"
#include "globals.hh"
#include "G4Run.hh"


MyRunAction::MyRunAction() : 
  theNumberOfEvents( 0 ), theSumOfTotalDepositedEnergy( 0.0 )
{}


MyRunAction::~MyRunAction() {}


void MyRunAction::BeginOfRunAction( const G4Run* aRun ) {
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  theNumberOfEvents = 0;
  theSumOfTotalDepositedEnergy = 0.0; 
}


void MyRunAction::EndOfRunAction( const G4Run* ) {
  double average = 0.0;
  if ( theNumberOfEvents > 0 ) {
    average =  theSumOfTotalDepositedEnergy / static_cast< double >( theNumberOfEvents );
  } 
  G4cout << G4endl << G4endl
         << " ###############  Final Statistics  ############### " << G4endl
         << "   Average deposited energy = " << average/GeV << " GeV" << G4endl
         << " ################################################## " << G4endl << G4endl;
}


void MyRunAction::updateEndOfEvent( const double eventTotalDepositedEnergy ) {
  theNumberOfEvents++;
  theSumOfTotalDepositedEnergy += eventTotalDepositedEnergy;
}

