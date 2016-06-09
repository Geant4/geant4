//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4WeightCutOffProcess.cc,v 1.7 2004/10/19 00:59:40 kurasige Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4WeightCutOffProcess.cc
//
// ----------------------------------------------------------------------

#include "G4WeightCutOffProcess.hh"
#include "G4VScorer.hh"
#include "G4GeometryCellStep.hh"
#include "G4ParallelGCellFinder.hh"
#include "G4MassGCellFinder.hh"
#include "G4TouchableHandle.hh"
#include "G4VIStore.hh"

G4WeightCutOffProcess::
G4WeightCutOffProcess(G4double wsurvival,
                      G4double wlimit,
                      G4double isource,
                      G4VIStore *istore,
                      const G4VGCellFinder &aGCellFinder,
                      const G4String &aName)
  : G4VProcess(aName), 
    fParticleChange(new G4ParticleChange),
    fWeightSurvival(wsurvival),
    fWeightLimit(wlimit),
    fSourceImportance(isource),
    fIStore(istore),
    fGCellFinder(aGCellFinder)
{
  if (!fParticleChange)
  {
    G4Exception("G4WeightCutOffProcess::G4WeightCutOffProcess()",
                "FatalError", FatalException,
                "Failed to allocate G4ParticleChange !");
  }
  G4VProcess::pParticleChange = fParticleChange;
}

G4WeightCutOffProcess::~G4WeightCutOffProcess()
{
  delete fParticleChange;
}

G4double G4WeightCutOffProcess::
PostStepGetPhysicalInteractionLength(const G4Track &,
                                     G4double, G4ForceCondition* condition)
{
  *condition = Forced;
  return kInfinity;
}
  
G4VParticleChange * 
G4WeightCutOffProcess::PostStepDoIt(const G4Track& aTrack, 
                                    const G4Step &aStep)
{
  fParticleChange->Initialize(aTrack);
  G4GeometryCell postCell = fGCellFinder.GetPostGeometryCell(aStep);
  G4double R = fSourceImportance;
  if (fIStore)
  {
    G4double i = fIStore->GetImportance(postCell);
    if (i>0)
    {
      R/=i;
    }
  }
  G4double w = aTrack.GetWeight();
  if (w<R*fWeightLimit)
  {
    G4double ws = fWeightSurvival*R;
    G4double p = w/(ws);
    if (G4UniformRand()<p)
    {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
    }
    else
    {
      fParticleChange->ProposeWeight(ws);
    }                  
  }
  return fParticleChange;
}

const G4String &G4WeightCutOffProcess::GetName() const
{
  return theProcessName;
}

G4double G4WeightCutOffProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&,
                                      G4double  ,
                                      G4double  ,
                                      G4double& ,
                                      G4GPILSelection*)
{
  return -1.0;
}

  
G4double G4WeightCutOffProcess::
AtRestGetPhysicalInteractionLength(const G4Track& ,
                                   G4ForceCondition*)
{
  return -1.0;
}
  
G4VParticleChange*
G4WeightCutOffProcess::AtRestDoIt(const G4Track&, const G4Step&)
{
  return 0;
}

G4VParticleChange*
G4WeightCutOffProcess::AlongStepDoIt(const G4Track&, const G4Step&)
{
  return 0;
}
