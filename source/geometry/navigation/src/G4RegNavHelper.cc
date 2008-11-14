#include "G4RegNavHelper.hh"

std::vector< std::pair<G4int,G4double> > G4RegNavHelper::theStepLengths;

void G4RegNavHelper::ClearStepLengths()
{
  theStepLengths.clear();
}
void G4RegNavHelper::AddStepLength(G4int copyNo, G4double slen )
{
  theStepLengths.push_back( std::pair<G4int,G4double>(copyNo,slen) );
}
