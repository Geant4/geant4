#include "G4IsoParticleChange.hh"

void G4IsoParticleChange::
Copy(G4VParticleChange * aResult)
{
  SetTrueStepLength(aResult->GetTrueStepLength());
  SetLocalEnergyDeposit(aResult->GetLocalEnergyDeposit());
  SetStatusChange(aResult->GetStatusChange());
  SetSteppingControl(aResult->GetSteppingControl());
  SetNumberOfSecondaries(aResult->GetNumberOfSecondaries());
  SetVerboseLevel(aResult->GetVerboseLevel());
  for(G4int i=0; i<aResult->GetNumberOfSecondaries(); i++)
  {
    AddSecondary(aResult->GetSecondary(i));
  }
}
