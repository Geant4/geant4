// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelStepCut.cc                               *
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R.Nartallo
// - First implementation of xray_telescope StepCut
// - Based on Chandra and XMM models
// 
//
// **********************************************************************

#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

#include "XrayTelStepCut.hh"

XrayTelStepCut::XrayTelStepCut(const G4String& aName)
  : G4VDiscreteProcess(aName),MaxChargedStep(DBL_MAX)
{
  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created "<< G4endl;
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

