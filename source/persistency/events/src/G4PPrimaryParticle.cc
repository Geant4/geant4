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
// $Id: G4PPrimaryParticle.cc,v 1.3 2001/07/11 10:02:16 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#include "G4PPrimaryParticle.hh"

#include "G4PrimaryParticle.hh"

G4PPrimaryParticle::G4PPrimaryParticle(G4PrimaryParticle* aParticle)
{
  PDGcode = aParticle->GetPDGcode();

  G4ParticleDefinition* ag4code  = aParticle->GetG4code();
  // if( ag4code != 0 )
  // { ...will look up G4PParticleDefinition and set the G4code...; }

  Px      = aParticle->GetPx(); 
  Py      = aParticle->GetPy(); 
  Pz      = aParticle->GetPz(); 

  G4PrimaryParticle* aNext = aParticle->GetNext();
  if( aNext != 0 )
  { nextParticle = new(ooThis()) G4PPrimaryParticle( aNext ); }

  G4PrimaryParticle* aDaughter = aParticle->GetDaughter();
  if( aDaughter != 0 )
  { daughterParticle = new(ooThis()) G4PPrimaryParticle( aDaughter ); }

  trackID = aParticle->GetTrackID();
  mass    = aParticle->GetMass();
  polX    = aParticle->GetPolX();
  polY    = aParticle->GetPolY();
  polZ    = aParticle->GetPolZ();
}

G4PPrimaryParticle::~G4PPrimaryParticle()
{
  if(nextParticle != 0)
  { HepDelete(nextParticle); }
  if(daughterParticle != 0)
  { HepDelete(daughterParticle); }
}

G4PrimaryParticle* G4PPrimaryParticle::MakeTransientObject()
{
  G4PrimaryParticle* aParticle = new G4PrimaryParticle(PDGcode, Px, Py, Pz);

  // aParticle->SetG4Code(...will lookup G4ParticleDefinition and set it...;)

  if(nextParticle != 0)
  { aParticle->SetNext( nextParticle->MakeTransientObject() ); }
  if(daughterParticle != 0)
  { aParticle->SetDaughter( daughterParticle->MakeTransientObject() ); }

  aParticle->SetTrackID( trackID );
  aParticle->SetMass( mass );
  aParticle->SetPolarization( polX, polY, polZ );

  return aParticle;
}

