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
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History:
//      17 August 2004  P. Gumplinger, T. MacPhail
// ------------------------------------------------------------
//
#include "G4DecayWithSpin.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4DecayTable.hh"
#include "G4MuonDecayChannelWithSpin.hh"

#include "G4Vector3D.hh"

#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"

#include "G4Transform3D.hh"

G4DecayWithSpin::G4DecayWithSpin(const G4String& processName):G4Decay(processName){}

G4DecayWithSpin::~G4DecayWithSpin(){}

G4VParticleChange* G4DecayWithSpin::DecayIt(const G4Track& aTrack, const G4Step& aStep)
{

// get particle
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();

// get parent_polarization
  G4ThreeVector parent_polarization = aParticle->GetPolarization();

  if(parent_polarization == G4ThreeVector(0,0,0))
  {
    // Generate random polarization direction

    G4double cost = 1. - 2.*G4UniformRand();
    G4double sint = std::sqrt((1.-cost)*(1.+cost));

    G4double phi = twopi*G4UniformRand();
    G4double sinp = std::sin(phi);
    G4double cosp = std::cos(phi);

    G4double px = sint*cosp;
    G4double py = sint*sinp;
    G4double pz = cost;

    parent_polarization.setX(px);
    parent_polarization.setY(py);
    parent_polarization.setZ(pz);

  }else{

    G4FieldManager* fieldMgr = aStep.GetTrack()->GetVolume()->
                                     GetLogicalVolume()->GetFieldManager();

    if (!fieldMgr) {
       G4TransportationManager *transportMgr =
                         G4TransportationManager::GetTransportationManager();
       G4PropagatorInField* fFieldPropagator = 
                                        transportMgr->GetPropagatorInField();
       if (fFieldPropagator) fieldMgr = 
                                  fFieldPropagator->GetCurrentFieldManager();
    }

    const G4Field* field = NULL;
    if(fieldMgr)field = fieldMgr->GetDetectorField();

    if (field && !(fieldMgr->DoesFieldChangeEnergy())) {

       G4double point[4];
       point[0] = (aStep.GetPostStepPoint()->GetPosition())[0];
       point[1] = (aStep.GetPostStepPoint()->GetPosition())[1];
       point[2] = (aStep.GetPostStepPoint()->GetPosition())[2];
       point[3] = aTrack.GetGlobalTime();

       G4double fieldValue[3];
       field -> GetFieldValue(point,fieldValue);

       G4ThreeVector B(fieldValue[0],fieldValue[1],fieldValue[2]);

       parent_polarization = Spin_Precession(aStep,B,fRemainderLifeTime);

    }
  }

// decay table
  G4DecayTable *decaytable = aParticleDef->GetDecayTable();

  if (decaytable) {
     G4MuonDecayChannelWithSpin *decaychannel;
     decaychannel = (G4MuonDecayChannelWithSpin*)decaytable->SelectADecayChannel();
     if (decaychannel) decaychannel->SetPolarization(parent_polarization);
  }

  G4ParticleChangeForDecay* pParticleChangeForDecay;

  pParticleChangeForDecay = (G4ParticleChangeForDecay*)G4Decay::DecayIt(aTrack,aStep);

  pParticleChangeForDecay->ProposePolarization(parent_polarization);

  return pParticleChangeForDecay;

}

G4ThreeVector G4DecayWithSpin::Spin_Precession( const G4Step& aStep,
                                        G4ThreeVector B, G4double deltatime )
{
  G4double Bnorm = std::sqrt(sqr(B[0])  + sqr(B[1]) +sqr(B[2]) );

  G4double q = aStep.GetTrack()->GetDefinition()->GetPDGCharge();
  G4double a = 1.165922e-3;
  G4double s_omega = 8.5062e+7*rad/(s*kilogauss);

  G4double omega = -(q*s_omega)*(1.+a) * Bnorm;

  G4double rotationangle = deltatime * omega;

  G4Transform3D SpinRotation = G4Rotate3D(rotationangle,B.unit());

  G4Vector3D Spin = aStep.GetTrack() -> GetPolarization();

  G4Vector3D newSpin = SpinRotation * Spin;

#ifdef G4VERBOSE

//  G4double normspin = std::sqrt(Spin*Spin);
//  G4double normnewspin = std::sqrt(newSpin*newSpin);
//  G4double cosalpha = Spin*newSpin/normspin/normnewspin;
//  G4double alpha = std::acos(cosalpha);

//  G4cout<< "AT REST::: PARAMETERS\n"
//        << "Initial spin  : " << Spin <<"\n"
//        << "Delta time    : " << deltatime <<"\n"
//        << "Rotation angle: " << rotationangle/rad <<"\n"
//        << "New spin      : " << newSpin <<"\n"
//        << "Checked norms : " << normspin <<" " << normnewspin <<" \n"
//        << G4endl;

#endif

  return newSpin;

}
