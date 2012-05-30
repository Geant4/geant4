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
// $Id: G4CollisionNN.cc,v 1.4 2010-03-12 15:45:18 gunter Exp $ //


#include "globals.hh"
#include "G4CollisionNN.hh"
#include "G4CollisionComposite.hh"
#include "G4VCollision.hh"
#include "G4CollisionVector.hh"
#include "G4KineticTrack.hh"
#include "G4VCrossSectionSource.hh"
#include "G4XNNTotal.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4CollisionNNElastic.hh"
#include "G4CollisionnpElastic.hh"
#include "G4CollisionNNToNDelta.hh"
#include "G4CollisionNNToDeltaDelta.hh"
#include "G4CollisionNNToNDeltastar.hh"
#include "G4CollisionNNToDeltaDeltastar.hh"
#include "G4CollisionNNToNNstar.hh"
#include "G4CollisionNNToDeltaNstar.hh"
#include "G4Pair.hh"

typedef GROUP8(G4CollisionnpElastic, G4CollisionNNElastic, 
              G4CollisionNNToNDelta, G4CollisionNNToDeltaDelta, 
	      G4CollisionNNToNDeltastar, G4CollisionNNToDeltaDeltastar,
	      G4CollisionNNToNNstar, G4CollisionNNToDeltaNstar) theChannels;

G4CollisionNN::G4CollisionNN()
{ 
  components=0;
  crossSectionSource = new G4XNNTotal();
  G4CollisionComposite::Register aR;
  G4ForEach<theChannels>::Apply(&aR, this);
}


G4CollisionNN::~G4CollisionNN()
{ 
  if (components) {
	delete components;
	components=0;
  }
  delete crossSectionSource;
  crossSectionSource = 0;
}


const std::vector<G4String>& G4CollisionNN::GetListOfColliders(G4int ) const
{
	  throw G4HadronicException(__FILE__, __LINE__, "G4CollisionNN::GetListOfColliders - Argument outside valid range"); 
	  return colliders1;
}


G4double G4CollisionNN::CrossSection(const G4KineticTrack& aTrk1, 
				    const G4KineticTrack& aTrk2) const
{
  G4double sigma = 0.;

  // nucleon-nucleon cross-sections made for on-shell particles.
  // here we take the kinetic energy as the quantity relevant
  // for calculating the scattering cross-sections for off-shell hadrons
  
  const G4VCrossSectionSource* xSource = GetCrossSectionSource();
  G4LorentzVector p1 = aTrk1.Get4Momentum();
  G4LorentzVector p2 = aTrk2.Get4Momentum();
  G4double t1 = p1.e()-aTrk1.GetActualMass();
  G4double t2 = p2.e()-aTrk2.GetActualMass();
  p1.setE(t1+aTrk1.GetDefinition()->GetPDGMass());
  p2.setE(t2+aTrk2.GetDefinition()->GetPDGMass());
  G4KineticTrack trk1(aTrk1);
  trk1.Set4Momentum(p1);
  G4KineticTrack trk2(aTrk2);
  trk2.Set4Momentum(p2);
  if( (p1+p2).mag()<aTrk1.GetDefinition()->GetPDGMass()+aTrk2.GetDefinition()->GetPDGMass())
  {
    return 0.;
  }

  if (xSource != 0)
    {
      // There is a cross section for this Collision
      sigma = xSource->CrossSection(trk1,trk2);
    }
  return sigma;
}

