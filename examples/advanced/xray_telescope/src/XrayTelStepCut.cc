//  XrayTelStepCut.cc

#include "XrayTelStepCut.hh"

#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

XrayTelStepCut::XrayTelStepCut(const G4String& aName)
  : G4VDiscreteProcess(aName),MaxChargedStep(DBL_MAX)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< endl;
   }
}

XrayTelStepCut::~XrayTelStepCut()
{
}

XrayTelStepCut::XrayTelStepCut(XrayTelStepCut& right)
    :G4VDiscreteProcess(right)
{}

void XrayTelStepCut::SetMaxStep(G4double step)
{
  MaxChargedStep = step ;
}

