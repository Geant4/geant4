#include "G4ImportanceTally.hh"

G4ImportanceTally::G4ImportanceTally(const G4String &tallyname) :
  fName(tallyname),
  fImportance(-1),
  fHasTallied(false)
{}


void G4ImportanceTally::Tally(const G4String &rawtallyname, 
			      G4Sigma &,
			      G4double importance){
  fImportance = importance;
  fHasTallied = true;
}
void G4ImportanceTally::Reset(){
  fImportance = -1;
  fHasTallied = false;
}
G4String G4ImportanceTally::GetName(){
  return fName;
}
G4double G4ImportanceTally::GetValue(){
  return fImportance;
}
G4bool G4ImportanceTally::HasTallied(){
  return fHasTallied;
}

