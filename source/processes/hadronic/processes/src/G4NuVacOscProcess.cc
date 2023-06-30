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
// Geant4 Hadron Elastic Scattering Process 
// 
// Created  from G4HadronElasticProcess
//  
// Modified:
//
// 5.4.23 V.Grichine - first implementation
//

#include <iostream>
#include <typeinfo>
#include "G4NuVacOscProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4Nucleus.hh"
#include "G4ProcessManager.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4ProductionCutsTable.hh"
#include "G4HadronicException.hh"
#include "G4HadronicInteraction.hh"
#include "G4VCrossSectionRatio.hh"
#include "G4VDiscreteProcess.hh"
#include "G4MuNeutrinoNucleusTotXsc.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "G4DynamicParticle.hh"
#include "G4StepPoint.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4SafetyHelper.hh"
#include "G4TransportationManager.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4NeutrinoTau.hh"
#include "G4PhysicsModelCatalog.hh"

///////////////////////////////////////////////////////////////////////////////


G4NuVacOscProcess::G4NuVacOscProcess(const G4String& eName, const G4String& pName)
  : G4HadronicProcess( pName, fHadronInelastic )
{
  fLowestEnergy = 1.*eV;
  fEnvelopeName = eName;
  theNuE       = G4NeutrinoE::NeutrinoE();
  theAntiNuE   = G4AntiNeutrinoE::AntiNeutrinoE();
  theNuMu      = G4NeutrinoMu::NeutrinoMu();
  theAntiNuMu  = G4AntiNeutrinoMu::AntiNeutrinoMu();
  theNuTau     = G4NeutrinoTau::NeutrinoTau();
  theAntiNuTau = G4AntiNeutrinoTau::AntiNeutrinoTau();

  InitParameters();
}

/////////////////////////////////////////////////////////
//
// Init the neutrino oscillation parameters

void G4NuVacOscProcess::InitParameters()
{
  if( fNormOrd ) // normal mass ordering
  {
    fSin2t12 = 0.31; 
    fSin2t23 = 0.558; 
    fSin2t13 = 0.02241;
    fDsm21   = 7.390e-5*CLHEP::eV*CLHEP::eV; 
    fDsm32   = 2.449e-3*CLHEP::eV*CLHEP::eV;
    fdcp     = CLHEP::degree * 222.; // 270.; // 90.; // 120.; //
  }
  else
  {
    fSin2t12 = 0.31; 
    fSin2t23 = 0.563; 
    fSin2t13 = 0.02261;
    fDsm21   = 7.3900e-5*CLHEP::eV*CLHEP::eV; 
    fDsm32   = -2.509e-3*CLHEP::eV*CLHEP::eV;
    fdcp     = CLHEP::degree * 285.; // 120. //
  }
  G4double c12(1.), s12(0.), c13(1.), s13(0.), c23(1.), s23(0.);

  s12 = std::sqrt( fSin2t12 );
  s23 = std::sqrt( fSin2t23 );
  s13 = std::sqrt( fSin2t13 );

  c12 = std::sqrt( 1. - fSin2t12 );
  c23 = std::sqrt( 1. - fSin2t23 );
  c13 = std::sqrt( 1. - fSin2t13 );

  G4complex expdcp = G4complex( std::cos(fdcp), std::sin(fdcp) ); // exp(i*deltaCP)

  G4complex u11, u12, u13, u21, u22, u23, u31, u32, u33;

  u11 =  c12*c13;                      u12 = c13*s12;                       u13 = s13*conj(expdcp);

  u21 = -s12*c23 - s13*s23*c12*expdcp; u22 = c12*c23 - s12*s23*s13*expdcp;  u23 = c13*s23;

  u31 = s12*s23 - s13*c12*c23*expdcp;  u32 = -c12*s23 - s12*s13*c23*expdcp; u33 = c13*c23;

  // fUdcp[3][3] = {  { u11, u12, u13 }, { u21, u22, u23 }, { u31, u32, u33 }  };

  // fUdcp[3][3] = {  u11, u12, u13,  u21, u22, u23,  u31, u32, u33  };
 
  fUdcp[0][0] = u11;  fUdcp[0][1] = u12;  fUdcp[0][2] = u13;
  fUdcp[1][0] = u21;  fUdcp[1][1] = u22;  fUdcp[1][2] = u23;
  fUdcp[2][0] = u31;  fUdcp[2][1] = u32;  fUdcp[2][2] = u33;

  G4double  m12, m13, m21, m23, m31, m32; //m11(0.), m22(0.),, m33(0.)

  m12 = -fDsm21; m13 = -fDsm21-fDsm32;
  m21 = -m12;    m23 = -fDsm32;
  m31 = -m13;    m32 = -m23;

  fDms[0][0] = fDms[1][1] = fDms[2][2] = 0.; // asymmetric
  fDms[0][1] = m12;         fDms[0][2] = m13;
  fDms[1][0] = m21;         fDms[1][2] = m23;
  fDms[2][0] = m31;         fDms[2][1] = m32;
}

////////////////////////////////////////////////////////////////////
//
// In long volumes, the mean free path can be reduced together with 
// the volume size along the neutrino trajectory

G4double G4NuVacOscProcess::
GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *)
{
  const G4String rName = 
    aTrack.GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetRegion()->GetName();
  G4double lambda(0.);
  G4double energy = aTrack.GetKineticEnergy();

  lambda = 0.4*CLHEP::hbarc*energy/( fDsm32 + fDsm21 );

  if( rName == fEnvelopeName && fNuNuclTotXscBias > 1.)    
  {
    lambda /= fNuNuclTotXscBias;
  }
  return lambda;
}

///////////////////////////////////////////////////

void G4NuVacOscProcess::ProcessDescription(std::ostream& outFile) const
{
  outFile << "G4NuVacOscProcess handles the oscillation of \n"
	  << "three flavor neutrinos on electrons by invoking the following  model(s) and \n"
	  << "mean pathe much smaller than the oscillation period.\n";
}

///////////////////////////////////////////////////////////////////////

G4VParticleChange* 
G4NuVacOscProcess::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  if( track.GetTrackStatus() != fAlive ) 
  { 
    return theTotalResult; 
  }
  theTotalResult->Clear();
  theTotalResult->Initialize(track);

  G4double weight = track.GetWeight();
  theTotalResult->ProposeWeight(weight);
  G4double kineticEnergy = track.GetKineticEnergy();

  if ( kineticEnergy <= fLowestEnergy )   return theTotalResult;

  const G4DynamicParticle*    dynParticle = track.GetDynamicParticle();
  const G4ParticleDefinition* part = dynParticle->GetDefinition();
  G4LorentzVector lv1 = dynParticle->Get4Momentum();
  
  G4int aa(0), bb(0); // neutrino flavors
  G4double ll = track.GetTrackLength(); // total track length
  const G4String rName =
    step.GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetRegion()->GetName();
  if(rName ==  fEnvelopeName && fNuNuclTotXscBias > 1.) ll *= fNuNuclTotXscBias;
  G4DynamicParticle* aLept =  nullptr;

  if( part == theAntiNuE  || 
      part == theAntiNuMu || 
      part == theAntiNuTau  ) fAnti = true;
  else                          fAnti = false;

  if(      part == theNuE   ||  part == theAntiNuE  )   aa = 0;
  else if( part == theNuMu  ||  part == theAntiNuMu  )  aa = 1;
  else                                                     aa = 2;

  bb = NuVacProbability( aa, kineticEnergy, ll); // oscillation engine

  if( bb == aa ) // no change
  {
    return theTotalResult;
  }
  else if( bb == 0 ) // new flavor (anti)neutrino - kill initial & add new
  {
    if( !fAnti ) aLept = new G4DynamicParticle( theNuE, lv1 );
    else         aLept = new G4DynamicParticle( theAntiNuE, lv1 );
  }
  else if( bb == 1 )
  {
    if( !fAnti ) aLept = new G4DynamicParticle( theNuMu, lv1 );
    else         aLept = new G4DynamicParticle( theAntiNuMu, lv1 );
  }
  else if( bb == 2 )
  {
    if( !fAnti ) aLept = new G4DynamicParticle( theNuTau, lv1 );
    else         aLept = new G4DynamicParticle( theAntiNuTau, lv1 );
  }
  theTotalResult->ProposeTrackStatus( fStopAndKill );
  theTotalResult->AddSecondary( aLept );

  return theTotalResult;
}

/////////////////////////////////////////////////////
//
// Oscillation probability aa->bb for neutrino energy Enu and its track distance Lnu

G4int G4NuVacOscProcess::NuVacProbability( G4int aa, G4double Enu, G4double Lnu )
{
  G4double probab(0.), probac(0.), probaa(0.), rr(0.), elCof(0.), delta[3][3];

  G4int bb(0), cc(0);

  if     ( aa == 0 ) { bb = 1; cc = 2; }
  else if( aa == 1 ) { bb = 0; cc = 2; }
  else if( aa == 2 ) { bb = 0; cc = 1; }

  elCof = 0.5*Lnu/Enu/CLHEP::hbarc;

  G4complex tmp(0.,0.), sum1(0.,0.),  sum2(0.,0.), expdel;

  for( G4int i = 0; i < 3; ++i )
  {
    for( G4int j = 0; j < 3; ++j ) delta[i][j] = fDms[i][j]*elCof;
  } 
  if( !fAnti )
  {
    for( G4int j = 0; j < 3; ++j )
    {
      for( G4int k = j+1; k < 3; ++k )
      {
        expdel = G4complex( std::cos( delta[k][j] ), -std::sin( delta[k][j] ) );

        tmp    = conj( fUdcp[bb][k] ) * fUdcp[aa][k] * fUdcp[bb][j] * conj( fUdcp[aa][j] );

        sum1  += tmp * std::sin( delta[k][j]*0.5 ) * std::sin( delta[k][j]*0.5 );
        sum2  += tmp * std::sin( delta[k][j] );
      }
    }
    probab = 2.*imag(sum2)  - 4.*real(sum1);

    sum1 = sum2 = G4complex( 0., 0. );

    for( G4int j = 0; j < 3; ++j )
    {
      for( G4int k = j+1; k < 3; ++k )
      {
        expdel = G4complex( std::cos( delta[k][j] ), -std::sin( delta[k][j] ) );

        tmp    = conj( fUdcp[cc][k] ) * fUdcp[aa][k] * fUdcp[cc][j] * conj( fUdcp[aa][j] );

        sum1  += tmp * std::sin( delta[k][j]*0.5 ) * std::sin( delta[k][j]*0.5 );
        sum2  += tmp * std::sin( delta[k][j] );
      }
    }
    probac = 2.*imag(sum2)  - 4.*real(sum1);
  }
  else // anti CP: exp(-i*delta)
  {
    for( G4int j = 0; j < 3; ++j )
    {
      for( G4int k = j+1; k < 3; ++k )
      {
        expdel = G4complex( std::cos( delta[k][j] ), -std::sin( delta[k][j] ) );

        tmp    = fUdcp[bb][k] * conj( fUdcp[aa][k] ) *conj( fUdcp[bb][j] ) * fUdcp[aa][j];

        sum1  += tmp * std::sin( delta[k][j]*0.5 ) * std::sin( delta[k][j]*0.5 );
        sum2  += tmp * std::sin( delta[k][j] );
      }
    }
    probab = 2.*imag(sum2)  - 4.*real(sum1);
    sum1 = sum2 = G4complex( 0., 0. );

    for( G4int j = 0; j < 3; ++j )
    {
      for( G4int k = j+1; k < 3; ++k )
      {
        expdel = G4complex( std::cos( delta[k][j] ), -std::sin( delta[k][j] ) );

        tmp    =  fUdcp[cc][k] * conj( fUdcp[aa][k] ) * conj( fUdcp[cc][j] ) * fUdcp[aa][j];

        sum1  += tmp * std::sin( delta[k][j]*0.5 ) * std::sin( delta[k][j]*0.5 );
        sum2  += tmp * std::sin( delta[k][j] );
      }
    }
    probac = 2.*imag(sum2)  - 4.*real(sum1);
  }
  probaa = 1. - probab - probac;

  if ( probaa < 0.) 
  {
    G4cout<<" sum neutrino disappearance > 1. "<<G4endl;

    rr = G4UniformRand()*( probab + probac );

    if( rr <= probab ) return bb;
    else return cc;
  }
  else
  {
    rr = G4UniformRand();

    if     ( rr <= probab )                          return bb;
    else if( rr >  probab && rr <= probab + probac ) return cc;
    else                                             return aa;
  }
}

//////////////////////////////////////////////////////////
