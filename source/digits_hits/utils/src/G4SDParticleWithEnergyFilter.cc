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
// $Id: G4SDParticleWithEnergyFilter.cc,v 1.1 2005/11/16 23:04:04 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// G4VSensitiveDetector
#include "G4SDParticleWithEnergyFilter.hh"
#include "G4Step.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

////////////////////////////////////////////////////////////////////////////////
// class description:
//
//  This is the class of a filter to be associated with a
// sensitive detector. 
//  This class filters steps by partilce definition and kinetic energy.
//
// Created: 2005-11-14  Tsukasa ASO.
// 
///////////////////////////////////////////////////////////////////////////////

G4SDParticleWithEnergyFilter::G4SDParticleWithEnergyFilter(G4String name,
							   G4double elow,
							   G4double ehigh)
  :G4VSDFilter(name)
{
  fParticleFilter = new G4SDParticleFilter(name);
  fKineticFilter  = new G4SDKineticEnergyFilter(name,elow,ehigh);
}

G4SDParticleWithEnergyFilter::~G4SDParticleWithEnergyFilter()
{ 
  delete fParticleFilter;
  delete fKineticFilter;
}

G4bool G4SDParticleWithEnergyFilter::Accept(const G4Step* aStep) const
{
  if ( ! fParticleFilter->Accept(aStep) )  return FALSE;
  if ( ! fKineticFilter->Accept(aStep)  )  return FALSE;
  return TRUE;
}

void G4SDParticleWithEnergyFilter::add(const G4String& particleName)
{
  fParticleFilter->add(particleName);
}

void G4SDParticleWithEnergyFilter::SetKineticEnergy(G4double elow, 
						    G4double ehigh)
{
  fKineticFilter->SetKineticEnergy(elow,ehigh);
}

void G4SDParticleWithEnergyFilter::show(){
  fParticleFilter->show();
  fKineticFilter->show();
}
