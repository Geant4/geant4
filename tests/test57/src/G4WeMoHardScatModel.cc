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
// $Id: G4WeMoHardScatModel.cc,v 1.2 2009-07-31 15:31:32 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4WeMoHardScatModel
//
// Author:        V. Grichine based on CoulombScatteringModel 
//
// Creation date: 31.07.2009
//
// Modifications:
//
//
// Class Description:
//
// -------------------------------------------------------------------
//
////////////////////////////////////////////////////////////////////////

#include "G4WeMoHardScatModel.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Proton.hh"

//////////////////////////////////////////////////////////////////////////

using namespace std;

G4WeMoHardScatModel::G4WeMoHardScatModel(const G4String& nam)
  : G4eWeMoHardScatModel(nam)
{}

///////////////////////////////////////////////////////////////////////////

G4WeMoHardScatModel::~G4WeMoHardScatModel()
{}

///////////////////////////////////////////////////////////////////////
//
//

G4double G4WeMoHardScatModel::ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition* p,
				G4double kinEnergy, 
				G4double Z, 
				G4double, 
				G4double cutEnergy,
				G4double)
{
  SetupParticle(p);
  G4double ekin = std::max(fLowEnergyLimit, kinEnergy);
  SetupKinematic(ekin, cutEnergy);

  // save lab system kinematics

  G4double xtkin = fTkin;
  G4double xmom2 = fMom2;
  G4double xinvb = fInvBeta2;

  // CM system

  fiz            = G4int(Z);
  G4double m2   = fNistManager->GetAtomicMassAmu(fiz)*amu_c2;
  G4double etot = fTkin + fMass;
  G4double ptot = sqrt(fMom2);

  G4double m12  = fMass*fMass;
  G4double momCM= ptot*m2/sqrt(m12 + m2*m2 + 2.0*etot*m2);

  fMom2 = momCM*momCM;
  fTkin = sqrt(fMom2 + m12) - fMass;

  //invbeta2 = 1.0 +  m12/mom2;
 
  G4double fm = m2/(fMass + m2);
  fInvBeta2 = 1.0 +  m12*fm*fm/fMom2;

  SetupTarget(Z, fTkin);

  G4double xsec = CrossSectionPerAtom();

  // restore Lab system kinematics
  fTkin = xtkin;
  fMom2 = xmom2;
  fInvBeta2 = xinvb;

  return xsec;
}

//////////////////////////////////////////////////////////////////////////////
//
//

void G4WeMoHardScatModel::SampleSecondaries(
			       std::vector<G4DynamicParticle*>* fvect,
			       const G4MaterialCutsCouple* couple,
			       const G4DynamicParticle* dp,
			       G4double cutEnergy, 
			       G4double)
{
  G4double kinEnergy = dp->GetKineticEnergy();
  if(kinEnergy <= DBL_MIN) return;
  DefineMaterial(couple);
  SetupParticle(dp->GetDefinition());
  G4double ekin = std::max(fLowEnergyLimit, kinEnergy);
  SetupKinematic(ekin, cutEnergy);

  // Choose nucleus

  fCurrentElement = SelectRandomAtom( couple, fParticle, ekin, feCut, fTkin);

  G4double Z  = fCurrentElement->GetZ();
  fiz         = G4int(Z);
  G4int ia    = SelectIsotopeNumber(fCurrentElement);
  G4double m2 = theParticleTable->GetIonTable()->GetNucleusMass(fiz, ia);

  // CM system
  G4double etot = fTkin + fMass;
  G4double ptot = sqrt(fMom2);

  G4double momCM = ptot*m2/sqrt(fMass*fMass + m2*m2 + 2.0*etot*m2);
  fMom2 = momCM*momCM;
  G4double m12 = fMass*fMass;
  G4double eCM = sqrt(fMom2 + m12);

  // a correction for heavy projectile

  G4double fm = m2/(fMass + m2);
  fInvBeta2 = 1.0 +  m12*fm*fm/fMom2;

  // sample scattering angle in CM system

  SetupTarget(Z, eCM - fMass);

  G4double cost = SampleCosineTheta();
  G4double z1   = 1.0 - cost;
  if(z1 < 0.0) return;
  G4double sint = sqrt(z1*(1.0 + cost));
  G4double phi  = twopi * G4UniformRand();

  // kinematics in the Lab system
  G4double bet  = ptot/(etot + m2);
  G4double gam  = 1.0/sqrt((1.0 - bet)*(1.0 + bet));
  G4double pzCM = momCM*cost;

  G4ThreeVector v1(momCM*cos(phi)*sint,momCM*sin(phi)*sint,gam*(pzCM + bet*eCM));
  G4ThreeVector dir = dp->GetMomentumDirection(); 
  G4ThreeVector newDirection = v1.unit();
  newDirection.rotateUz(dir);   
  fParticleChange->ProposeMomentumDirection(newDirection);   

  //  G4double elab = gam*(eCM + bet*pzCM);

  G4double Ecm  = sqrt(fMass*fMass + m2*m2 + 2.0*etot*m2);
  G4double elab = etot - m2*(ptot/Ecm)*(ptot/Ecm)*(1.-cost) ;

 
  ekin = elab - fMass;
  if(ekin < 0.0) ekin = 0.0;
  fParticleChange->SetProposedKineticEnergy(ekin);

  // recoil
  G4double erec = kinEnergy - ekin;

  if(erec > fRecoilThreshold) 
  {
    G4ParticleDefinition* ion = theParticleTable->FindIon(fiz, ia, 0, fiz);
    G4double plab = sqrt(ekin*(ekin + 2.0*fMass));
    G4ThreeVector p2 = (ptot*dir - plab*newDirection).unit();
    G4DynamicParticle* newdp  = new G4DynamicParticle(ion, p2, erec);
    fvect->push_back(newdp);
  } 
  else if(erec > 0.0) 
  {
    fParticleChange->ProposeLocalEnergyDeposit(erec);
    fParticleChange->ProposeNonIonizingEnergyDeposit(erec);
  }
}

//
//
///////////////////////////////////////////////////////////////////////////////////


