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
// $Id: G4VElasticCollision.cc,v 1.3 2006-06-29 20:41:53 gunter Exp $ //

#include "globals.hh"
#include "G4VElasticCollision.hh"
#include "G4KineticTrack.hh"
#include "G4VCrossSectionSource.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4XNNElastic.hh"
#include "G4AngularDistribution.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4KineticTrackVector.hh"
#include "G4AngularDistributionNP.hh"   //  np scattering
#include "G4AngularDistributionPP.hh"   //  nn and pp scattering
#include <typeinfo>

G4VElasticCollision::G4VElasticCollision()
{ 
}


G4VElasticCollision::~G4VElasticCollision()
{ }


G4KineticTrackVector* G4VElasticCollision::FinalState(const G4KineticTrack& trk1,
						      const G4KineticTrack& trk2) const
{
  const G4VAngularDistribution* angDistribution;

  angDistribution = GetAngularDistribution();


  G4LorentzVector pCM=trk1.Get4Momentum() + trk2.Get4Momentum();

  G4LorentzRotation toLabFrame(pCM.boostVector());
  G4LorentzVector Ptmp=toLabFrame.inverse() * trk1.Get4Momentum();  //trk1 in CMS
  G4LorentzRotation toZ;
  toZ.rotateZ(-Ptmp.phi());
  toZ.rotateY(-Ptmp.theta());
  toLabFrame *= toZ.inverse();

  G4double S = pCM.mag2();
  G4double m10 = trk1.GetDefinition()->GetPDGMass();
  G4double m20 = trk2.GetDefinition()->GetPDGMass();
  if(S-(m10+m20)*(m10+m20) < 0) return new G4KineticTrackVector;

  G4double m_1 = trk1.GetActualMass();
  G4double m_2 = trk2.GetActualMass();
  
  // Angles of outgoing particles
  G4double cosTheta = angDistribution->CosTheta(S,m_1,m_2);

   if ( (trk1.GetDefinition() == G4Proton::Proton() || trk1.GetDefinition() == G4Neutron::Neutron() )
      &&(trk2.GetDefinition() == G4Proton::Proton() || trk2.GetDefinition() == G4Neutron::Neutron() ) )
   {
      if ( trk1.GetDefinition() == trk2.GetDefinition() )
      {
	  if ( trk1.GetDefinition() == G4Proton::Proton() )
	  {
//	      G4cout << "scatterangle pp " << cosTheta
//	             << " " << typeid(*angDistribution).name() << G4endl;
	  } else {
//	      G4cout << "scatterangle nn " << cosTheta
//	             << " " << typeid(*angDistribution).name() << G4endl;
	  }
      } else {
//	  G4cout << "scatterangle pn " << cosTheta
//	             << " " << typeid(*angDistribution).name() << G4endl;
      }
   } else {
//      G4cout << "scatterangle other " << cosTheta
//	             << " " << typeid(*angDistribution).name() << G4endl;
   }

  G4double phi = angDistribution->Phi();
  G4double Theta = std::acos(cosTheta);

  // Unit vector of three-momentum
  G4ThreeVector pFinal1(std::sin(Theta)*std::cos(phi), std::sin(Theta)*std::sin(phi), cosTheta);
  // Three momentum in cm system
  G4double pInCM = std::sqrt((S-(m10+m20)*(m10+m20))*(S-(m10-m20)*(m10-m20))/(4.*S));
  pFinal1 = pFinal1 * pInCM;
  G4ThreeVector pFinal2 = -pFinal1;

  G4double eFinal1 = std::sqrt(pFinal1.mag2() + m10*m10);
  G4double eFinal2 = std::sqrt(pFinal2.mag2() + m20*m20);

  G4LorentzVector p4Final1(pFinal1, eFinal1);
  G4LorentzVector p4Final2(pFinal2, eFinal2);

  // Lorentz transformation
  p4Final1 *= toLabFrame;
  p4Final2 *= toLabFrame;

  // Final tracks are copies of incoming ones, with modified 4-momenta
  G4KineticTrack* final1 = new G4KineticTrack(trk1);
  final1->Set4Momentum(p4Final1);
  G4KineticTrack* final2 = new G4KineticTrack(trk2);
  final2->Set4Momentum(p4Final2);

  G4KineticTrackVector* finalTracks = new G4KineticTrackVector;
  finalTracks->push_back(final1);
  finalTracks->push_back(final2);

  return finalTracks;
}
