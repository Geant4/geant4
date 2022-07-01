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
// $Id: G4NuMuNucleusCcModel.cc 91806 2015-08-06 12:20:45Z gcosmo $
//
// Geant4 Header : G4NuMuNucleusCcModel
//
// Author : V.Grichine 12.2.19
//  

#include <iostream>
#include <fstream>
#include <sstream>

#include "G4NuMuNucleusCcModel.hh"
// #include "G4NuMuNuclCcDistrKR.hh" 

// #include "G4NuMuResQX.hh" 

#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
// #include "G4Threading.hh"

// #include "G4Integrator.hh"
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"
/*
#include "G4CascadeInterface.hh"
// #include "G4BinaryCascade.hh"
#include "G4TheoFSGenerator.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4FTFModel.hh"
// #include "G4BinaryCascade.hh"
#include "G4HadFinalState.hh"
#include "G4HadSecondary.hh"
#include "G4HadronicInteractionRegistry.hh"
// #include "G4INCLXXInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSMFragmentation.hh"
#include "G4QGSParticipants.hh"
*/
#include "G4KineticTrack.hh"
#include "G4DecayKineticTracks.hh"
#include "G4KineticTrackVector.hh"
#include "G4Fragment.hh"
#include "G4NucleiProperties.hh"
#include "G4ReactionProductVector.hh"

#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"


#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Nucleus.hh"
#include "G4LorentzVector.hh"

using namespace std;
using namespace CLHEP;

#ifdef G4MULTITHREADED
    G4Mutex G4NuMuNucleusCcModel::numuNucleusModel = G4MUTEX_INITIALIZER;
#endif     


G4NuMuNucleusCcModel::G4NuMuNucleusCcModel(const G4String& name) 
  : G4NeutrinoNucleusModel(name)
{
  fData = fMaster = false;
  InitialiseModel();  
}


G4NuMuNucleusCcModel::~G4NuMuNucleusCcModel()
{}


void G4NuMuNucleusCcModel::ModelDescription(std::ostream& outFile) const
{

    outFile << "G4NuMuNucleusCcModel is a neutrino-nucleus (charge current)  scattering\n"
            << "model which uses the standard model \n"
            << "transfer parameterization.  The model is fully relativistic\n";

}

/////////////////////////////////////////////////////////
//
// Read data from G4PARTICLEXSDATA (locally PARTICLEXSDATA)

void G4NuMuNucleusCcModel::InitialiseModel()
{
  G4String pName  = "nu_mu";
  // G4String pName  = "anti_nu_mu";
  
  G4int nSize(0), i(0), j(0), k(0);

  if(!fData)
  { 
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&numuNucleusModel);
    if(!fData)
    { 
#endif     
      fMaster = true;
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&numuNucleusModel);
#endif
  }
  
  if(fMaster)
  {  
    const char* path = G4FindDataDir("G4PARTICLEXSDATA");
    std::ostringstream ost1, ost2, ost3, ost4;
    ost1 << path << "/" << "neutrino" << "/" << pName << "/xarraycckr";

    std::ifstream filein1( ost1.str().c_str() );

    // filein.open("$PARTICLEXSDATA/");

    filein1>>nSize;

    for( k = 0; k < fNbin; ++k )
    {
      for( i = 0; i <= fNbin; ++i )
      {
        filein1 >> fNuMuXarrayKR[k][i];
        // G4cout<< fNuMuXarrayKR[k][i] << "  ";
      }
    }
    // G4cout<<G4endl<<G4endl;

    ost2 << path << "/" << "neutrino" << "/" << pName << "/xdistrcckr";
    std::ifstream  filein2( ost2.str().c_str() );

    filein2>>nSize;

    for( k = 0; k < fNbin; ++k )
    {
      for( i = 0; i < fNbin; ++i )
      {
        filein2 >> fNuMuXdistrKR[k][i];
        // G4cout<< fNuMuXdistrKR[k][i] << "  ";
      }
    }
    // G4cout<<G4endl<<G4endl;

    ost3 << path << "/" << "neutrino" << "/" << pName << "/q2arraycckr";
    std::ifstream  filein3( ost3.str().c_str() );

    filein3>>nSize;

    for( k = 0; k < fNbin; ++k )
    {
      for( i = 0; i <= fNbin; ++i )
      {
        for( j = 0; j <= fNbin; ++j )
        {
          filein3 >> fNuMuQarrayKR[k][i][j];
          // G4cout<< fNuMuQarrayKR[k][i][j] << "  ";
        }
      }
    }
    // G4cout<<G4endl<<G4endl;

    ost4 << path << "/" << "neutrino" << "/" << pName << "/q2distrcckr";
    std::ifstream  filein4( ost4.str().c_str() );

    filein4>>nSize;

    for( k = 0; k < fNbin; ++k )
    {
      for( i = 0; i <= fNbin; ++i )
      {
        for( j = 0; j < fNbin; ++j )
        {
          filein4 >> fNuMuQdistrKR[k][i][j];
          // G4cout<< fNuMuQdistrKR[k][i][j] << "  ";
        }
      }
    }
    fData = true;
  }
}

/////////////////////////////////////////////////////////

G4bool G4NuMuNucleusCcModel::IsApplicable(const G4HadProjectile & aPart, 
					        G4Nucleus & )
{
  G4bool result  = false;
  G4String pName = aPart.GetDefinition()->GetParticleName();
  G4double energy = aPart.GetTotalEnergy();
  
  if(  pName == "nu_mu"  // || pName == "anti_nu_mu" )  
        &&
        energy > fMinNuEnergy                                )
  {
    result = true;
  }

  return result;
}

/////////////////////////////////////////// ClusterDecay ////////////////////////////////////////////////////////////
//
//

G4HadFinalState* G4NuMuNucleusCcModel::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();
  fProton = f2p2h = fBreak = false;
  fCascade = fString  = false;
  fLVh = fLVl = fLVt = fLVcpi = G4LorentzVector(0.,0.,0.,0.);

  const G4HadProjectile* aParticle = &aTrack;
  G4double energy = aParticle->GetTotalEnergy();

  G4String pName  = aParticle->GetDefinition()->GetParticleName();

  if( energy < fMinNuEnergy ) 
  {
    theParticleChange.SetEnergyChange(energy);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

  SampleLVkr( aTrack, targetNucleus);

  if( fBreak == true || fEmu < fMu ) // ~5*10^-6
  {
    // G4cout<<"ni, ";
    theParticleChange.SetEnergyChange(energy);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

  // LVs of initial state

  G4LorentzVector lvp1 = aParticle->Get4Momentum();
  G4LorentzVector lvt1( 0., 0., 0., fM1 );
  G4double mPip = G4ParticleTable::GetParticleTable()->FindParticle(211)->GetPDGMass();

  // 1-pi by fQtransfer && nu-energy
  G4LorentzVector lvpip1( 0., 0., 0., mPip );
  G4LorentzVector lvsum, lv2, lvX;
  G4ThreeVector eP;
  G4double cost(1.), sint(0.), phi(0.), muMom(0.), massX2(0.), massX(0.), massR(0.), eCut(0.);
  G4DynamicParticle* aLept = nullptr; // lepton lv

  G4int Z = targetNucleus.GetZ_asInt();
  G4int A = targetNucleus.GetA_asInt();
  G4double  mTarg = targetNucleus.AtomicMass(A,Z);
  G4int pdgP(0), qB(0);
  // G4double mSum = G4ParticleTable::GetParticleTable()->FindParticle(2212)->GetPDGMass() + mPip;

  G4int iPi     = GetOnePionIndex(energy);
  G4double p1pi = GetNuMuOnePionProb( iPi, energy);

  if( p1pi > G4UniformRand()  && fCosTheta > 0.9  ) // && fQtransfer < 0.95*GeV ) // mu- & coherent pion + nucleus
  {
    // lvsum = lvp1 + lvpip1;
    lvsum = lvp1 + lvt1;
    // cost = fCosThetaPi;
    cost = fCosTheta;
    sint = std::sqrt( (1.0 - cost)*(1.0 + cost) );
    phi  = G4UniformRand()*CLHEP::twopi;
    eP   = G4ThreeVector( sint*std::cos(phi), sint*std::sin(phi), cost );

    // muMom = sqrt(fEmuPi*fEmuPi-fMu*fMu);
    muMom = sqrt(fEmu*fEmu-fMu*fMu);

    eP *= muMom;

    // lv2 = G4LorentzVector( eP, fEmuPi );
    // lv2 = G4LorentzVector( eP, fEmu );
    lv2 = fLVl;

    // lvX = lvsum - lv2;
    lvX = fLVh;
    massX2 = lvX.m2();
    massX = lvX.m();
    massR = fLVt.m();
    
    if ( massX2 <= 0. ) // vmg: very rarely ~ (1-4)e-6 due to big Q2/x, to be improved
    {
      fCascade = true;
      theParticleChange.SetEnergyChange(energy);
      theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
      return &theParticleChange;
    }
    fW2 = massX2;

    if(  pName == "nu_mu" )         aLept = new G4DynamicParticle( theMuonMinus, lv2 );  
    // else if( pName == "anti_nu_mu") aLept = new G4DynamicParticle( theMuonPlus,  lv2 );
    else
    {
      theParticleChange.SetEnergyChange(energy);
      theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
      return &theParticleChange;
    }
    if( pName == "nu_mu" ) pdgP =  211;
    // else                   pdgP = -211;
    // eCut = fMpi + 0.5*(fMpi*fMpi-massX2)/mTarg; // massX -> fMpi

    if( A > 1 )
    {
      eCut = (fMpi + mTarg)*(fMpi + mTarg) - (massX + massR)*(massX + massR);
      eCut /= 2.*massR;
      eCut += massX;
    }
    else  eCut = fM1 + fMpi;

    if ( lvX.e() > eCut ) // && sqrt( GetW2() ) < 1.4*GeV ) // 
    {
      CoherentPion( lvX, pdgP, targetNucleus);
    }
    else
    {
      fCascade = true;
      theParticleChange.SetEnergyChange(energy);
      theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
      return &theParticleChange;
    } 
    theParticleChange.AddSecondary( aLept, fSecID );

    return &theParticleChange;
  }
  else // lepton part in lab
  { 
    lvsum = lvp1 + lvt1;
    cost = fCosTheta;
    sint = std::sqrt( (1.0 - cost)*(1.0 + cost) );
    phi  = G4UniformRand()*CLHEP::twopi;
    eP   = G4ThreeVector( sint*std::cos(phi), sint*std::sin(phi), cost );

    muMom = sqrt(fEmu*fEmu-fMu*fMu);

    eP *= muMom;

    lv2 = G4LorentzVector( eP, fEmu );
    lv2 = fLVl;
    lvX = lvsum - lv2;
    lvX = fLVh;
    massX2 = lvX.m2();

    if ( massX2 <= 0. ) // vmg: very rarely ~ (1-4)e-6 due to big Q2/x, to be improved
    {
      fCascade = true;
      theParticleChange.SetEnergyChange(energy);
      theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
      return &theParticleChange;
    }
    fW2 = massX2;

    if(  pName == "nu_mu" )         aLept = new G4DynamicParticle( theMuonMinus, lv2 );  
    // else if( pName == "anti_nu_mu") aLept = new G4DynamicParticle( theMuonPlus,  lv2 );
    else
    {
      theParticleChange.SetEnergyChange(energy);
      theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
      return &theParticleChange;
    }
    theParticleChange.AddSecondary( aLept, fSecID );
  }

  // hadron part

  fRecoil  = nullptr;
  
  if( A == 1 )
  {
    if( pName == "nu_mu" ) qB = 2;
    // else                   qB = 0;

    // if( G4UniformRand() > 0.1 ) //  > 0.9999 ) // > 0.0001 ) //
    {
      ClusterDecay( lvX, qB );
    }
    return &theParticleChange;
  }
    /*
    // else
    {
      if( pName == "nu_mu" ) pdgP =  211;
      else                   pdgP = -211;


      if ( fQtransfer < 0.95*GeV ) // < 0.35*GeV ) //
      {
	if( lvX.m() > mSum ) CoherentPion( lvX, pdgP, targetNucleus);
      }
    }
    return &theParticleChange;
  }
  */
  G4Nucleus recoil;
  G4double rM(0.), ratio = G4double(Z)/G4double(A);

  if( ratio > G4UniformRand() ) // proton is excited
  {
    fProton = true;
    recoil = G4Nucleus(A-1,Z-1);
    fRecoil = &recoil;
    rM = recoil.AtomicMass(A-1,Z-1);

    if( pName == "nu_mu" ) // (++) state -> p + pi+
    { 
      fMt = G4ParticleTable::GetParticleTable()->FindParticle(2212)->GetPDGMass()
          + G4ParticleTable::GetParticleTable()->FindParticle(211)->GetPDGMass();
    }
    else // (0) state -> p + pi-, n + pi0
    {
      // fMt = G4ParticleTable::GetParticleTable()->FindParticle(2212)->GetPDGMass()
      //     + G4ParticleTable::GetParticleTable()->FindParticle(-211)->GetPDGMass();
    } 
  }
  else // excited neutron
  {
    fProton = false;
    recoil = G4Nucleus(A-1,Z);
    fRecoil = &recoil;
    rM = recoil.AtomicMass(A-1,Z);

    if( pName == "nu_mu" ) // (+) state -> n + pi+
    {      
      fMt = G4ParticleTable::GetParticleTable()->FindParticle(2112)->GetPDGMass()
          + G4ParticleTable::GetParticleTable()->FindParticle(211)->GetPDGMass();
    }
    else // (-) state -> n + pi-, // n + pi0
    {
      // fMt = G4ParticleTable::GetParticleTable()->FindParticle(2112)->GetPDGMass()
      //     + G4ParticleTable::GetParticleTable()->FindParticle(-211)->GetPDGMass();
    } 
  }
  // G4int       index = GetEnergyIndex(energy);
  G4int nepdg = aParticle->GetDefinition()->GetPDGEncoding();

  G4double qeTotRat; // = GetNuMuQeTotRat(index, energy);
  qeTotRat = CalculateQEratioA( Z, A, energy, nepdg);

  G4ThreeVector dX = (lvX.vect()).unit();
  G4double eX   = lvX.e();  // excited nucleon
  G4double mX   = sqrt(massX2);
  // G4double pX   = sqrt( eX*eX - mX*mX );
  // G4double sumE = eX + rM;

  if( qeTotRat > G4UniformRand() || mX <= fMt ) // || eX <= 1232.*MeV) // QE
  {  
    fString = false;

    if( fProton ) 
    {  
      fPDGencoding = 2212;
      fMr =  proton_mass_c2;
      recoil = G4Nucleus(A-1,Z-1);
      fRecoil = &recoil;
      rM = recoil.AtomicMass(A-1,Z-1);
    } 
    else // if( pName == "anti_nu_mu" ) 
    {  
      fPDGencoding = 2112;
      fMr =   G4ParticleTable::GetParticleTable()->
	FindParticle(fPDGencoding)->GetPDGMass(); // 939.5654133*MeV;
      recoil = G4Nucleus(A-1,Z);
      fRecoil = &recoil;
      rM = recoil.AtomicMass(A-1,Z);
    } 
    // sumE = eX + rM;   
    G4double eTh = fMr + 0.5*(fMr*fMr - mX*mX)/rM;

    if( eX <= eTh ) // vmg, very rarely out of kinematics
    {
      fString = true;
      theParticleChange.SetEnergyChange(energy);
      theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
      return &theParticleChange;
    }
    // FinalBarion( fLVh, 0, fPDGencoding ); // p(n)+deexcited recoil
    FinalBarion( lvX, 0, fPDGencoding ); // p(n)+deexcited recoil
  }
  else // if ( eX < 9500000.*GeV ) // <  25.*GeV) // < 95.*GeV ) // < 2.5*GeV ) //cluster decay
  {  
    if     (  fProton && pName == "nu_mu" )      qB =  2;
    // else if(  fProton && pName == "anti_nu_mu" ) qB =  0;
    else if( !fProton && pName == "nu_mu" )      qB =  1;
    // else if( !fProton && pName == "anti_nu_mu" ) qB = -1;


      ClusterDecay( lvX, qB );
  }
  return &theParticleChange;
}


/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////
//
// sample x, then Q2

void G4NuMuNucleusCcModel::SampleLVkr(const G4HadProjectile & aTrack, G4Nucleus& targetNucleus)
{
  fBreak = false;
  G4int A = targetNucleus.GetA_asInt(), iTer(0), iTerMax(100); 
  G4int Z = targetNucleus.GetZ_asInt(); 
  G4double e3(0.), pMu2(0.), pX2(0.), nMom(0.), rM(0.), hM(0.), tM = targetNucleus.AtomicMass(A,Z);
  G4double Ex(0.), ei(0.), nm2(0.);
  G4double cost(1.), sint(0.), phi(0.), muMom(0.); 
  G4ThreeVector eP, bst;
  const G4HadProjectile* aParticle = &aTrack;
  G4LorentzVector lvp1 = aParticle->Get4Momentum();

  if( A == 1 ) // hydrogen, no Fermi motion ???
  {
    fNuEnergy = aParticle->GetTotalEnergy();
    iTer = 0;

    do
    {
      fXsample = SampleXkr(fNuEnergy);
      fQtransfer = SampleQkr(fNuEnergy, fXsample);
      fQ2 = fQtransfer*fQtransfer;

     if( fXsample > 0. )
      {
        fW2 = fM1*fM1 - fQ2 + fQ2/fXsample; // sample excited hadron mass
        fEmu = fNuEnergy - fQ2/2./fM1/fXsample;
      }
      else
      {
        fW2 = fM1*fM1;
        fEmu = fNuEnergy;
      }
      e3 = fNuEnergy + fM1 - fEmu;

      if( e3 < sqrt(fW2) )  G4cout<<"energyX = "<<e3/GeV<<", fW = "<<sqrt(fW2)/GeV<<G4endl;
    
      pMu2 = fEmu*fEmu - fMu*fMu;

      if(pMu2 < 0.) { fBreak = true; return; }

      pX2  = e3*e3 - fW2;

      fCosTheta  = fNuEnergy*fNuEnergy  + pMu2 - pX2;
      fCosTheta /= 2.*fNuEnergy*sqrt(pMu2);
      iTer++;
    }
    while( ( abs(fCosTheta) > 1. || fEmu < fMu ) && iTer < iTerMax );

    if( iTer >= iTerMax ) { fBreak = true; return; }

    if( abs(fCosTheta) > 1.) // vmg: due to big Q2/x values. To be improved ...
    { 
      G4cout<<"H2: fCosTheta = "<<fCosTheta<<", fEmu = "<<fEmu<<G4endl;
      // fCosTheta = -1. + 2.*G4UniformRand(); 
      if(fCosTheta < -1.) fCosTheta = -1.;
      if(fCosTheta >  1.) fCosTheta =  1.;
    }
    // LVs

    G4LorentzVector lvt1  = G4LorentzVector( 0., 0., 0., fM1 );
    G4LorentzVector lvsum = lvp1 + lvt1;

    cost = fCosTheta;
    sint = std::sqrt( (1.0 - cost)*(1.0 + cost) );
    phi  = G4UniformRand()*CLHEP::twopi;
    eP   = G4ThreeVector( sint*std::cos(phi), sint*std::sin(phi), cost );
    muMom = sqrt(fEmu*fEmu-fMu*fMu);
    eP *= muMom;
    fLVl = G4LorentzVector( eP, fEmu );

    fLVh = lvsum - fLVl;
    fLVt = G4LorentzVector( 0., 0., 0., 0. ); // no recoil
  }
  else // Fermi motion, Q2 in nucleon rest frame
  {
    G4Nucleus recoil1( A-1, Z );
    rM = recoil1.AtomicMass(A-1,Z);   
    do
    {
      // nMom = NucleonMomentumBR( targetNucleus ); // BR
      nMom = GgSampleNM( targetNucleus ); // Gg
      Ex = GetEx(A-1, fProton);
      ei = tM - sqrt( (rM + Ex)*(rM + Ex) + nMom*nMom );
      //   ei = 0.5*( tM - s2M - 2*eX );
    
      nm2 = ei*ei - nMom*nMom;
      iTer++;
    }
    while( nm2 < 0. && iTer < iTerMax ); 

    if( iTer >= iTerMax ) { fBreak = true; return; }
    
    G4ThreeVector nMomDir = nMom*G4RandomDirection();

    if( !f2p2h || A < 3 ) // 1p1h
    {
      // hM = tM - rM;

      fLVt = G4LorentzVector( -nMomDir, sqrt( (rM + Ex)*(rM + Ex) + nMom*nMom ) ); // rM ); //
      fLVh = G4LorentzVector(  nMomDir, ei ); // hM); //
    }
    else // 2p2h
    {
      G4Nucleus recoil(A-2,Z-1);
      rM = recoil.AtomicMass(A-2,Z-1)+sqrt(nMom*nMom+fM1*fM1);
      hM = tM - rM;

      fLVt = G4LorentzVector( nMomDir, sqrt( rM*rM+nMom*nMom ) );
      fLVh = G4LorentzVector(-nMomDir, sqrt( hM*hM+nMom*nMom )  ); 
    }
    // G4cout<<hM<<", ";
    // bst = fLVh.boostVector();

    // lvp1.boost(-bst); // -> nucleon rest system, where Q2 transfer is ???

    fNuEnergy  = lvp1.e();
    // G4double mN = fLVh.m(); // better mN = fM1 !? vmg
    iTer = 0;

    do // no FM!?, 5.4.20 vmg
    {
      fXsample = SampleXkr(fNuEnergy);
      fQtransfer = SampleQkr(fNuEnergy, fXsample);
      fQ2 = fQtransfer*fQtransfer;

      // G4double mR = mN + fM1*(A-1.)*std::exp(-2.0*fQtransfer/mN); // recoil mass in+el

      if( fXsample > 0. )
      {
        fW2 = fM1*fM1 - fQ2 + fQ2/fXsample; // sample excited hadron mass

        // fW2 = mN*mN - fQ2 + fQ2/fXsample; // sample excited hadron mass
        // fEmu = fNuEnergy - fQ2/2./mR/fXsample; // fM1->mN

        fEmu = fNuEnergy - fQ2/2./fM1/fXsample; // fM1->mN
      }
      else
      {
        // fW2 = mN*mN;

        fW2 = fM1*fM1; 
        fEmu = fNuEnergy;
      }
      // if(fEmu < 0.) G4cout<<"fEmu = "<<fEmu<<" hM = "<<hM<<G4endl;
      // e3 = fNuEnergy + mR - fEmu;

      e3 = fNuEnergy + fM1 - fEmu;

      // if( e3 < sqrt(fW2) )  G4cout<<"energyX = "<<e3/GeV<<", fW = "<<sqrt(fW2)/GeV<<G4endl;
    
      pMu2 = fEmu*fEmu - fMu*fMu;
      pX2  = e3*e3 - fW2;

      if(pMu2 < 0.) { fBreak = true; return; }

      fCosTheta  = fNuEnergy*fNuEnergy  + pMu2 - pX2;
      fCosTheta /= 2.*fNuEnergy*sqrt(pMu2);
      iTer++;
    }
    while( ( abs(fCosTheta) > 1. || fEmu < fMu ) && iTer < iTerMax );

    if( iTer >= iTerMax ) { fBreak = true; return; }

    if( abs(fCosTheta) > 1.) // vmg: due to big Q2/x values. To be improved ...
    { 
      G4cout<<"FM: fCosTheta = "<<fCosTheta<<", fEmu = "<<fEmu<<G4endl;
      // fCosTheta = -1. + 2.*G4UniformRand(); 
      if( fCosTheta < -1.) fCosTheta = -1.;
      if( fCosTheta >  1.) fCosTheta =  1.;
    }
    // LVs
    // G4LorentzVector lvt1  = G4LorentzVector( 0., 0., 0., mN ); // fM1 );

    G4LorentzVector lvt1  = G4LorentzVector( 0., 0., 0., fM1 ); // fM1 );
    G4LorentzVector lvsum = lvp1 + lvt1;

    cost = fCosTheta;
    sint = std::sqrt( (1.0 - cost)*(1.0 + cost) );
    phi  = G4UniformRand()*CLHEP::twopi;
    eP   = G4ThreeVector( sint*std::cos(phi), sint*std::sin(phi), cost );
    muMom = sqrt(fEmu*fEmu-fMu*fMu);
    eP *= muMom;
    fLVl = G4LorentzVector( eP, fEmu );
    fLVh = lvsum - fLVl;

    // if( fLVh.e() < mN || fLVh.m2() < 0.) { fBreak = true; return; }

    if( fLVh.e() < fM1 || fLVh.m2() < 0.) { fBreak = true; return; }

    // back to lab system

    // fLVl.boost(bst);
    // fLVh.boost(bst);
  }
  //G4cout<<iTer<<", "<<fBreak<<"; ";
}

//
//
///////////////////////////
