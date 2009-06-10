//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4VLowEnergyDiscretePhotonProcess.cc,v 1.6 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------
//
// File name:     G4VLowEnergyDiscretePhotonProcess.cc
//
// Author:        Capra Riccardo
//
// Creation date: May 2005
//
// History:
// -----------
// 02 May 2005  R. Capra         1st implementation
//
//----------------------------------------------------------------

  

#include "G4VLowEnergyDiscretePhotonProcess.hh"

#include "G4String.hh"
#include "G4CrossSectionHandler.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4Gamma.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh" // G4UniformRand

G4VLowEnergyDiscretePhotonProcess :: G4VLowEnergyDiscretePhotonProcess(const G4String& processName, 
								       const G4String& aCrossSectionFileName, 
								       const G4String& aScatterFileName, 
								       G4VDataSetAlgorithm* aScatterInterpolation, 
								       G4double aLowEnergyLimit, 
								       G4double aHighEnergyLimit)
  :
  G4VLowEnergyTestableDiscreteProcess(processName),
  lowEnergyLimit(aLowEnergyLimit),
  highEnergyLimit(aHighEnergyLimit),
  crossSectionFileName(aCrossSectionFileName),
  meanFreePathTable(0)
{
  crossSectionHandler = new G4CrossSectionHandler();
  scatterFunctionData = new G4CompositeEMDataSet(aScatterInterpolation, 1., 1.);
  scatterFunctionData->LoadData(aScatterFileName);
 
  if (verboseLevel > 0)
    {
      G4cout << GetProcessName() << " is created " << G4endl
	     << "Energy range: "
	     << lowEnergyLimit / keV << " keV - "
	     << highEnergyLimit / GeV << " GeV"
	     << G4endl;
    }
}



G4VLowEnergyDiscretePhotonProcess::~G4VLowEnergyDiscretePhotonProcess(void)
{
  if (meanFreePathTable)
    delete meanFreePathTable;
 
  delete crossSectionHandler;
  delete scatterFunctionData;
}





G4bool G4VLowEnergyDiscretePhotonProcess::IsApplicable(const G4ParticleDefinition& particleDefinition)
{
  return (&particleDefinition)==G4Gamma::Gamma(); 
}



void G4VLowEnergyDiscretePhotonProcess::BuildPhysicsTable(const G4ParticleDefinition&  /*photon*/)
{
  crossSectionHandler->Clear();
  crossSectionHandler->LoadData(crossSectionFileName);

  if (meanFreePathTable)
    delete meanFreePathTable; 
  meanFreePathTable=crossSectionHandler->BuildMeanFreePathForMaterials();
}





G4double G4VLowEnergyDiscretePhotonProcess::GetMeanFreePath(const G4Track& aTrack, G4double /*previousStepSize*/, G4ForceCondition*  /*condition*/)
{
  G4double photonEnergy;
  photonEnergy = aTrack.GetDynamicParticle()->GetKineticEnergy();

  if (photonEnergy < lowEnergyLimit)
    return DBL_MAX;
 
  if (photonEnergy > highEnergyLimit)
    photonEnergy=highEnergyLimit;
 
  size_t materialIndex;
  materialIndex = aTrack.GetMaterialCutsCouple()->GetIndex(); 

  return meanFreePathTable->FindValue(photonEnergy, materialIndex);
}





G4ThreeVector G4VLowEnergyDiscretePhotonProcess::GetPhotonPolarization(const G4DynamicParticle&  photon)
{
  G4ThreeVector photonMomentumDirection;
  G4ThreeVector photonPolarization;

  photonPolarization = photon.GetPolarization(); 
  photonMomentumDirection = photon.GetMomentumDirection();

  if ((!photonPolarization.isOrthogonal(photonMomentumDirection, 1e-6)) || photonPolarization.mag()==0.)
    {
      // if |photonPolarization|==0. or |photonPolarization * photonDirection0| > 1e-6 * |photonPolarization ^ photonDirection0|
      // then polarization is choosen randomly.
  
      G4ThreeVector e1(photonMomentumDirection.orthogonal().unit());
      G4ThreeVector e2(photonMomentumDirection.cross(e1).unit());
  
      G4double angle(G4UniformRand() * twopi);
  
      e1*=std::cos(angle);
      e2*=std::sin(angle);
  
      photonPolarization=e1+e2;
    }
  else if (photonPolarization.howOrthogonal(photonMomentumDirection) != 0.)
    {
      // if |photonPolarization * photonDirection0| != 0.
      // then polarization is made orthonormal;
  
      photonPolarization=photonPolarization.perpPart(photonMomentumDirection);
    }
 
  return photonPolarization.unit();
}
