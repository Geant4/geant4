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
//
// -- Bogus -- BaBar Object-Oriented Geant-based Unified Simulation
//
// NTSTLooperDeath
//
// Description:
//   This is a simple GEANT4 process that destroys any particle below
//   the specified total kinetic energy that reverse direction (loops 180
//   degress) in the x/y plane.
//
//   Based on BgsChargedLowEnergyDeath
//
// Author List:
//   David Williams
//
// Modification History:
//
//-----------------------------------------------------------------------------

#include "NTSTLooperDeath.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"

//
// Constructor
//
NTSTLooperDeath::NTSTLooperDeath( G4double theMinMomentum,
				  const char* name,
				  G4ProcessType type )
  : G4VProcess( name, type ), minMomentum( theMinMomentum )
{;}

NTSTLooperDeath::~NTSTLooperDeath(){;}
//
// PostStepGetPhysicalInteractionLength
//
G4double NTSTLooperDeath::PostStepGetPhysicalInteractionLength( 
						const G4Track& track,
					        G4double  , // previousStepSize,
					        G4ForceCondition* condition ) 
{
  const G4DynamicParticle *particle = track.GetDynamicParticle();
	
  *condition = NotForced;

  //
  // We don't touch any particle above the cut momentum
  //
  if (particle->GetTotalMomentum() > minMomentum) return DBL_MAX;
  
  //
  // Nor do we touch any particle with small transverse component
  // to their momentum. Here the cutoff is somewhat arbitrary.
  //
  G4double vperp = track.GetMomentumDirection().perp();
  if (vperp < 0.1) return DBL_MAX;
  
  //
  // How far in azimuthal angle do we need to go before the
  // particle loops back on itself?
  //
  G4ThreeVector initialDir(track.GetVertexMomentumDirection());
  G4ThreeVector dx = track.GetPosition() - track.GetVertexPosition();
  G4double dxPerp = dx.perp();
  if (dxPerp < 1E-6) return DBL_MAX;
  
  G4double dot = ( initialDir.x()*dx.x() + initialDir.y()*dx.y() )/dxPerp;
  
  if (dot < 0) return 0;		// Already done so
  
  G4double phi = pi - std::acos(dot);
  
  if (phi < 0) return 0;

  //
  // What is the radius of curvature?
  //
  // Only use the z component of the field to calculate this.
  //
  G4FieldManager *fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  
  if (!fieldManager->DoesFieldExist()) return DBL_MAX;
  
  //
  // Assume the field is a magnetic field (i.e. returns
  // a vector of three doubles of unit Telsa).
  //
  // There is no good way I know to confirm this, so it is purely
  // a matter of faith. *** BEWARE ***
  //
  G4MagneticField *field = (G4MagneticField *)fieldManager->GetDetectorField();
  G4double b[3];
  G4ThreeVector pos(track.GetPosition());
  G4double posv[3] = { pos.x(), pos.y(), pos.z() };
  field->GetFieldValue( posv, b );
  
  //
  // No field? Forget it!
  //
  if (std::fabs(b[2]) < 0.00001) return DBL_MAX;
  
  //
  // Calculate radius of curvature, the usual way. 
  // Note G4 default units: mm, Telsa, MeV
  //
  // G4 suggestion: the constant below should be added to
  // the geant4 list.
  //
  G4double radius = std::fabs(track.GetMomentum().perp()/299.79251/b[2]);
  
  //
  // Convert this to a distance
  //
  return radius*phi/vperp;
}


//
// PostStepDoit
//
G4VParticleChange *
NTSTLooperDeath::PostStepDoIt( const G4Track &track, const G4Step & ) // step ) 
{
  pParticleChange->Initialize(track);

  //
  // This is a tough one. What should happen when we kill off a looper?
  // For now: just deposit remaining energy (including mass).
  //
  const G4DynamicParticle *particle = track.GetDynamicParticle();
  G4double energyDeposited = particle->GetTotalEnergy();

  pParticleChange->ProposeTrackStatus( fStopAndKill );
  pParticleChange->SetNumberOfSecondaries( 0 );
  pParticleChange->ProposeLocalEnergyDeposit( energyDeposited );
  ClearNumberOfInteractionLengthLeft();

  return pParticleChange;
}






