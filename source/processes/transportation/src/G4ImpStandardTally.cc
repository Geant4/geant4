#include "G4ImpStandardTally.hh"
#include "G4StandardTally.hh"
#include "G4Sigma.hh"

G4ImpStandardTally::G4ImpStandardTally(const G4String &tallyname,
				 const G4String &rawtallyname,
				 const G4String &SigmaSpec) : 
  fStandardTally(new G4StandardTally(tallyname, rawtallyname, SigmaSpec))
{}

void G4ImpStandardTally::Reset(){
  fStandardTally->Reset();
}

void G4ImpStandardTally::Tally(const G4String &rawtallyname, 
			       G4Sigma &sigma, G4double importance){
  fStandardTally->Tally(rawtallyname, sigma);
  fValue = fStandardTally->GetValue() * importance;
}

G4String G4ImpStandardTally::GetName(){
  return fStandardTally->GetName();
}
G4double G4ImpStandardTally::GetValue(){
  return fValue;
}
G4bool G4ImpStandardTally::HasTallied(){
  return fStandardTally->HasTallied();
}

