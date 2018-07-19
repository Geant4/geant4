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
// $Id: G4NeutronElectronElModel.cc 91806 2015-08-06 12:20:45Z gcosmo $
//
// Geant4 Header : G4NeutronElectronElModel
//
//  16.5.17: V.Grichine
//  

#include "G4NeutronElectronElModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4Integrator.hh"
#include "G4Electron.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsFreeVector.hh"


using namespace std;
using namespace CLHEP;

G4NeutronElectronElModel::G4NeutronElectronElModel(const G4String& name) 
  : G4HadronElastic(name)
{
 // neutron magneton squared

  fM   = neutron_mass_c2; // neutron mass
  fM2  = fM*fM;
  fme  = electron_mass_c2;
  fme2 = fme*fme;
  fMv2 = 0.7056*GeV*GeV;

  SetMinEnergy( 0.001*GeV );
  SetMaxEnergy( 10.*TeV );
  SetLowestEnergyLimit(1.e-6*eV);  

  theElectron = G4Electron::Electron();
  // PDG2016: sin^2 theta Weinberg

  fEnergyBin = 200;
  fMinEnergy = 1.*MeV;
  fMaxEnergy = 10000.*GeV;
  fEnergyVector = new G4PhysicsLogVector(fMinEnergy, fMaxEnergy, fEnergyBin);

  fAngleBin = 500;
  fAngleTable = 0;

  fCutEnergy = 0.; // default value

  Initialise();
}

////////////////////////////////////////////////

G4NeutronElectronElModel::~G4NeutronElectronElModel()
{
  if( fEnergyVector ) 
  {
    delete fEnergyVector;
    fEnergyVector = 0;
  }
  if( fAngleTable )
  {
    fAngleTable->clearAndDestroy();
    delete fAngleTable;
    fAngleTable = nullptr;
  }
}

/////////////////////////////////////////

void G4NeutronElectronElModel::ModelDescription(std::ostream& outFile) const
{

    outFile << "G4NeutronElectronElModel is a neutrino-electron (neutral current) elastic scattering\n"
            << "model which uses the standard model \n"
            << "transfer parameterization.  The model is fully relativistic\n";

}

/////////////////////////////////////////////////////////

G4bool G4NeutronElectronElModel::IsApplicable(const G4HadProjectile & aTrack, 
  			      G4Nucleus & targetNucleus)
{
  G4bool result  = false;
  G4String pName = aTrack.GetDefinition()->GetParticleName();
  // G4double minEnergy = 0.; 
  G4double energy = aTrack.GetTotalEnergy();

  if( fCutEnergy > 0. ) // min detected recoil electron energy
  {
    // minEnergy = 0.5*(fCutEnergy+sqrt(fCutEnergy*(fCutEnergy+2.*electron_mass_c2)));
  }
  if( pName == "neutron"   &&
      energy >= fMinEnergy  && energy <= fMaxEnergy   )                            
  {
    result = true;
  }
  G4int Z = targetNucleus.GetZ_asInt();
        Z *= 1;

  return result;
}

////////////////////////////////////////////////////

void  G4NeutronElectronElModel::Initialise()
{
  G4double result = 0., sum, Tkin, dt, t1, t2;
  G4int iTkin, jTransfer;
  G4Integrator<G4NeutronElectronElModel, G4double(G4NeutronElectronElModel::*)(G4double)> integral;

  fAngleTable = new G4PhysicsTable(fEnergyBin);

  for( iTkin = 0; iTkin < fEnergyBin; iTkin++)
  {
    Tkin  = fEnergyVector->GetLowEdgeEnergy(iTkin);
    fAm      = CalculateAm(Tkin);
    dt = 1./fAngleBin;

    G4PhysicsFreeVector* vectorT = new G4PhysicsFreeVector(fAngleBin);

    sum = 0.;

    for( jTransfer = 0; jTransfer < fAngleBin; jTransfer++)
    {
      t1 = dt*jTransfer;
      t2 = t1 + dt;

      result = integral.Legendre96( this, &G4NeutronElectronElModel::XscIntegrand, t1, t2 );

      sum += result;
      // G4cout<<sum<<", ";
      vectorT->PutValue(jTransfer, t1, sum);
    }
    // G4cout<<G4endl;   
    fAngleTable->insertAt(iTkin,vectorT);
  }
  return;
}

//////////////////////////////////////////////////////
//
// sample recoil electron energy in lab frame

G4double G4NeutronElectronElModel::SampleSin2HalfTheta(G4double Tkin)
{
  G4double result = 0., position; 
  G4int iTkin, iTransfer;

  for( iTkin = 0; iTkin < fEnergyBin; iTkin++)
  {
    if( Tkin < fEnergyVector->GetLowEdgeEnergy(iTkin) ) break;
  }  
  if ( iTkin >= fEnergyBin ) iTkin = fEnergyBin-1;   // Tkin is more then theMaxEnergy
  if ( iTkin < 0 )           iTkin = 0; // against negative index, Tkin < theMinEnergy

    position = (*(*fAngleTable)(iTkin))(fAngleBin-1)*G4UniformRand();

    // G4cout<<"position = "<<position<<G4endl;

    for( iTransfer = 0; iTransfer < fAngleBin; iTransfer++)
    {
      if( position <= (*(*fAngleTable)(iTkin))(iTransfer) ) break;
    }
    if (iTransfer >= fAngleBin-1) iTransfer = fAngleBin-1;

    // G4cout<<"iTransfer = "<<iTransfer<<G4endl;

    result = GetTransfer(iTkin, iTransfer, position);

    // G4cout<<"t = "<<t<<G4endl;
  

  return result;
}

/////////////////////////////////////////////////

G4double 
G4NeutronElectronElModel:: GetTransfer( G4int iTkin, G4int iTransfer, G4double position )
{
  G4double x1, x2, y1, y2, randTransfer, delta, mean, epsilon = 1.e-6;

  if( iTransfer == 0 ||  iTransfer == fAngleBin-1 )
  {
    randTransfer = (*fAngleTable)(iTkin)->GetLowEdgeEnergy(iTransfer);
    // iTransfer++;
  }
  else
  {
    if ( iTransfer >= G4int((*fAngleTable)(iTkin)->GetVectorLength()) )
    {
      iTransfer = (*fAngleTable)(iTkin)->GetVectorLength() - 1;
    }
    y1 = (*(*fAngleTable)(iTkin))(iTransfer-1);
    y2 = (*(*fAngleTable)(iTkin))(iTransfer);

    x1 = (*fAngleTable)(iTkin)->GetLowEdgeEnergy(iTransfer-1);
    x2 = (*fAngleTable)(iTkin)->GetLowEdgeEnergy(iTransfer);

    delta = y2 - y1;
    mean  = y2 + y1;

    if ( x1 == x2 ) randTransfer = x2;
    else
    {
      // if ( y1 == y2 ) 

      if ( delta < epsilon*mean ) 
      {
        randTransfer = x1 + ( x2 - x1 )*G4UniformRand();
      }
      else 
      {
        randTransfer = x1 + ( position - y1 )*( x2 - x1 )/delta; // ( y2 - y1 );
      }
    }
  }
  return randTransfer;
}

//////////////////////////////////////////////////////////////
//
// Rosenbluth relation (ultra-relativistic!) in the neutron rest frame, 
// x = sin^2(theta/2), theta is the electron scattering angle
// Magnetic form factor in the dipole approximation.

G4double G4NeutronElectronElModel::XscIntegrand(G4double x)
{
  G4double result = 1., q2, znq2, znf, znf2, znf4;

  znq2 = 1. + 2.*fee*x/fM;

  q2 = 4.*fee2*x/znq2;

  znf  = 1 + q2/fMv2;
  znf2 = znf*znf;
  znf4 = znf2*znf2;

  result /= ( x + fAm )*znq2*znq2*znf4; 

  result *= ( 1 - x )/( 1 + q2/4./fM2 ) + 2.*x;

  return result;
}

////////////////////////////////////////////////
//
//

G4HadFinalState* G4NeutronElectronElModel::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();

  const G4HadProjectile* aParticle = &aTrack;
  G4double Tkin = aParticle->GetKineticEnergy();
  fAm = CalculateAm( Tkin);
  //   G4double En = aParticle->GetTotalEnergy();

  if( Tkin <= LowestEnergyLimit() ) 
  {
    theParticleChange.SetEnergyChange(Tkin);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }
  // sample e-scattering angle and make final state in lab frame

  G4double sin2ht = SampleSin2HalfTheta( Tkin); // in n-rrest frame

  // G4cout<<"sin2ht = "<<sin2ht<<G4endl;

  G4double eTkin = fee; // fM;

  eTkin /= 1.+2.*fee*sin2ht/fM; // fme/En + 2*sin2ht;

  eTkin -= fme;

  // G4cout<<"eTkin = "<<eTkin<<G4endl;

  if( eTkin > fCutEnergy )
  {
    G4double ePlab = sqrt( eTkin*(eTkin + 2.*fme) );

    // G4cout<<"ePlab = "<<ePlab<<G4endl;

    G4double cost = 1. - 2*sin2ht;

    if( cost >  1. ) cost = 1.;
    if( cost < -1. ) cost = -1.;

    G4double sint = std::sqrt( (1.0 - cost)*(1.0 + cost) );
    G4double phi  = G4UniformRand()*CLHEP::twopi;

    G4ThreeVector eP( sint*std::cos(phi), sint*std::sin(phi), cost );
    eP *= ePlab;
    G4LorentzVector lvt2( eP, eTkin + electron_mass_c2 ); // recoil e- in n-rest frame

    G4LorentzVector lvp1 = aParticle->Get4Momentum();
    G4LorentzVector lvt1(0.,0.,0.,electron_mass_c2);
    G4LorentzVector lvsum = lvp1+lvt1;

    G4ThreeVector bst = lvp1.boostVector();
    lvt2.boost(bst);

    // G4cout<<"lvt2 = "<<lvt2<<G4endl;

    G4DynamicParticle * aSec = new G4DynamicParticle( theElectron, lvt2 );
    theParticleChange.AddSecondary( aSec );

    G4LorentzVector lvp2 = lvsum-lvt2;

    // G4cout<<"lvp2 = "<<lvp2<<G4endl;

    G4double Tkin2 = lvp2.e()-aParticle->GetDefinition()->GetPDGMass();
    theParticleChange.SetEnergyChange(Tkin2);
    theParticleChange.SetMomentumChange(lvp2.vect().unit());
  }
  else if( eTkin > 0.0 ) 
  {
    theParticleChange.SetLocalEnergyDeposit( eTkin );
    Tkin -= eTkin;

    if( Tkin > 0. )
    {
      theParticleChange.SetEnergyChange( Tkin );
      theParticleChange.SetMomentumChange( aTrack.Get4Momentum().vect().unit() );
    }
  }
  else 
  {
    theParticleChange.SetEnergyChange( Tkin );
    theParticleChange.SetMomentumChange( aTrack.Get4Momentum().vect().unit() );
  }
  G4int Z = targetNucleus.GetZ_asInt();
        Z *= 1;
 
  return &theParticleChange;
}

//
//
///////////////////////////
