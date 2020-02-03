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
// $Id: G4NuMuNucleusNcModel.cc 91806 2015-08-06 12:20:45Z gcosmo $
//
// Geant4 Header : G4NuMuNucleusNcModel
//
// Author : V.Grichine 12.2.19
//  

#include "G4NuMuNucleusNcModel.hh"
#include "G4NeutrinoNucleusModel.hh" 

// #include "G4NuMuResQX.hh" 

#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"

// #include "G4Integrator.hh"
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"
#include "G4KineticTrack.hh"
#include "G4DecayKineticTracks.hh"
#include "G4KineticTrackVector.hh"
#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"


#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4Nucleus.hh"
#include "G4LorentzVector.hh"

using namespace std;
using namespace CLHEP;

const G4int G4NuMuNucleusNcModel::fResNumber = 6;

const G4double G4NuMuNucleusNcModel::fResMass[6] = // [fResNumber] = 
  {2190., 1920., 1700., 1600., 1440., 1232. };

const G4int G4NuMuNucleusNcModel::fClustNumber = 4;

const G4double G4NuMuNucleusNcModel::fMesMass[4] = {1260., 980., 770., 139.57};
const G4int    G4NuMuNucleusNcModel::fMesPDG[4]  = {20213, 9000211, 213, 211};

// const G4double G4NuMuNucleusNcModel::fBarMass[4] = {1905., 1600., 1232., 939.57};
// const G4int    G4NuMuNucleusNcModel::fBarPDG[4]  = {2226, 32224, 2224, 2212};

const G4double G4NuMuNucleusNcModel::fBarMass[4] = {1700., 1600., 1232., 939.57};
const G4int    G4NuMuNucleusNcModel::fBarPDG[4]  = {12224, 32224, 2224, 2212};

const G4double  G4NuMuNucleusNcModel::fNuMuEnergyLogVector[50] = {
115.603, 133.424, 153.991, 177.729, 205.126, 236.746, 273.24, 315.361, 363.973, 420.08, 484.836, 559.573, 645.832, 
745.387, 860.289, 992.903, 1145.96, 1322.61, 1526.49, 1761.8, 2033.38, 2346.83, 2708.59, 3126.12, 3608.02, 4164.19, 
4806.1, 5546.97, 6402.04, 7388.91, 8527.92, 9842.5, 11359.7, 13110.8, 15131.9, 17464.5, 20156.6, 23263.8, 26849.9, 
30988.8, 35765.7, 41279, 47642.2, 54986.3, 63462.4, 73245.2, 84536, 97567.2, 112607, 129966 };

G4double G4NuMuNucleusNcModel::fNuMuXarrayKR[50][51] = {{1.0}};
G4double G4NuMuNucleusNcModel::fNuMuXdistrKR[50][50] = {{1.0}};
G4double G4NuMuNucleusNcModel::fNuMuQarrayKR[50][51][51] = {{{1.0}}};
G4double G4NuMuNucleusNcModel::fNuMuQdistrKR[50][51][50] = {{{1.0}}};

#ifdef G4MULTITHREADED
    G4Mutex G4NuMuNucleusNcModel::numuNucleusModel = G4MUTEX_INITIALIZER;
#endif     


G4NuMuNucleusNcModel::G4NuMuNucleusNcModel(const G4String& name) 
  : G4NeutrinoNucleusModel(name)
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( 100.*TeV );
  SetMinEnergy(1.e-6*eV); 

  fMnumu = 0.; 
  fData = fMaster = false;
  InitialiseModel();  
     
}


G4NuMuNucleusNcModel::~G4NuMuNucleusNcModel()
{}


void G4NuMuNucleusNcModel::ModelDescription(std::ostream& outFile) const
{

    outFile << "G4NuMuNucleusNcModel is a neutrino-nucleus (neutral current) scattering\n"
            << "model which uses the standard model \n"
            << "transfer parameterization.  The model is fully relativistic\n";

}

/////////////////////////////////////////////////////////
//
// Read data from G4PARTICLEXSDATA (locally PARTICLEXSDATA)

void G4NuMuNucleusNcModel::InitialiseModel()
{
  G4String pName  = "nu_mu";
  
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
    char* path = getenv("G4PARTICLEXSDATA");
    std::ostringstream ost1, ost2, ost3, ost4;
    ost1 << path << "/" << "neutrino" << "/" << pName << "/xarraynckr";

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

    ost2 << path << "/" << "neutrino" << "/" << pName << "/xdistrnckr";
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

    ost3 << path << "/" << "neutrino" << "/" << pName << "/q2arraynckr";
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

    ost4 << path << "/" << "neutrino" << "/" << pName << "/q2distrnckr";
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

G4bool G4NuMuNucleusNcModel::IsApplicable(const G4HadProjectile & aPart, 
					       G4Nucleus & targetNucleus)
{
  G4bool result  = false;
  G4String pName = aPart.GetDefinition()->GetParticleName();
  G4double energy = aPart.GetTotalEnergy();
  
  if(  pName == "nu_mu" // || pName == "anti_nu_mu"   ) 
        &&
        energy > fMinNuEnergy                                )
  {
    result = true;
  }
  G4int Z = targetNucleus.GetZ_asInt();
        Z *= 1;

  return result;
}

/////////////////////////////////////////// ClusterDecay ////////////////////////////////////////////////////////////
//
//

G4HadFinalState* G4NuMuNucleusNcModel::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();
  fProton = f2p2h = fBreak = false;
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

  if( fBreak == true || fEmu < fMnumu ) // ~5*10^-6
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
  G4double cost(1.), sint(0.), phi(0.), muMom(0.), massX2(0.);
  G4DynamicParticle* aLept = nullptr; // lepton lv

  G4int Z = targetNucleus.GetZ_asInt();
  G4int A = targetNucleus.GetA_asInt();
  G4double  mTarg = targetNucleus.AtomicMass(A,Z);
  G4int pdgP(0), qB(0);
  // G4double mSum = G4ParticleTable::GetParticleTable()->FindParticle(2212)->GetPDGMass() + mPip;

  G4int iPi     = GetOnePionIndex(energy);
  G4double p1pi = GetNuMuOnePionProb( iPi, energy);

  if( p1pi > G4UniformRand()  ) // && fQtransfer < 0.95*GeV ) // mu- & coherent pion + nucleus
  {
    // lvsum = lvp1 + lvpip1;
    lvsum = lvp1 + lvt1;
    // cost = fCosThetaPi;
    cost = fCosTheta;
    sint = std::sqrt( (1.0 - cost)*(1.0 + cost) );
    phi  = G4UniformRand()*CLHEP::twopi;
    eP   = G4ThreeVector( sint*std::cos(phi), sint*std::sin(phi), cost );

    // muMom = sqrt(fEmuPi*fEmuPi-fMnumu*fMnumu);
    muMom = sqrt(fEmu*fEmu-fMnumu*fMnumu);

    eP *= muMom;

    // lv2 = G4LorentzVector( eP, fEmuPi );
    lv2 = G4LorentzVector( eP, fEmu );
    lv2 = fLVl;

    lvX = lvsum - lv2;
    lvX = fLVh;
    massX2 = lvX.m2();

    if ( massX2 <= 0. ) // vmg: very rarely ~ (1-4)e-6 due to big Q2/x, to be improved
    {
      theParticleChange.SetEnergyChange(energy);
      theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
      return &theParticleChange;
    }
    fW2 = massX2;

    if(  pName == "nu_mu" )         aLept = new G4DynamicParticle( theNuMu, lv2 );  
    else if( pName == "anti_nu_mu") aLept = new G4DynamicParticle( theANuMu,  lv2 );
    else
    {
      theParticleChange.SetEnergyChange(energy);
      theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
      return &theParticleChange;
    } 
 
    pdgP = 111;

    G4double eCut = fMpi + 0.5*(fMpi*fMpi - massX2)/mTarg; // massX -> fMpi

    if ( lvX.e() > eCut ) // && sqrt( GetW2() ) < 1.4*GeV ) // 
    {
      CoherentPion( lvX, pdgP, targetNucleus);
    }
    else
    {
      theParticleChange.SetEnergyChange(energy);
      theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
      return &theParticleChange;
    } 
    theParticleChange.AddSecondary( aLept );

    return &theParticleChange;
  }
  else // lepton part in lab
  { 
    lvsum = lvp1 + lvt1;
    cost = fCosTheta;
    sint = std::sqrt( (1.0 - cost)*(1.0 + cost) );
    phi  = G4UniformRand()*CLHEP::twopi;
    eP   = G4ThreeVector( sint*std::cos(phi), sint*std::sin(phi), cost );

    muMom = sqrt(fEmu*fEmu-fMnumu*fMnumu);

    eP *= muMom;

    lv2 = G4LorentzVector( eP, fEmu );

    lvX = lvsum - lv2;

    massX2 = lvX.m2();

    if ( massX2 <= 0. ) // vmg: very rarely ~ (1-4)e-6 due to big Q2/x, to be improved
    {
      theParticleChange.SetEnergyChange(energy);
      theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
      return &theParticleChange;
    }
    fW2 = massX2;

    if(  pName == "nu_mu" )         aLept = new G4DynamicParticle( theNuMu, lv2 );  
    else if( pName == "anti_nu_mu") aLept = new G4DynamicParticle( theANuMu,  lv2 );
    
    theParticleChange.AddSecondary( aLept );
  }

  // hadron part

  fRecoil  = nullptr;
  fCascade = false;
  fString  = false;
  
  if( A == 1 )
  {
    qB = 1;

    // if( G4UniformRand() > 0.1 ) //  > 0.9999 ) // > 0.0001 ) //
    {
      ClusterDecay( lvX, qB );
    }
    return &theParticleChange;
  }
  G4Nucleus recoil;
  G4double rM(0.), ratio = G4double(Z)/G4double(A);

  if( ratio > G4UniformRand() ) // proton is excited
  {
    fProton = true;
    recoil = G4Nucleus(A-1,Z-1);
    fRecoil = &recoil;
    rM = recoil.AtomicMass(A-1,Z-1);

    fMt = G4ParticleTable::GetParticleTable()->FindParticle(2212)->GetPDGMass()
          + G4ParticleTable::GetParticleTable()->FindParticle(111)->GetPDGMass();
  }
  else // excited neutron
  {
    fProton = false;
    recoil = G4Nucleus(A-1,Z);
    fRecoil = &recoil;
    rM = recoil.AtomicMass(A-1,Z);

    fMt = G4ParticleTable::GetParticleTable()->FindParticle(2112)->GetPDGMass()
          + G4ParticleTable::GetParticleTable()->FindParticle(111)->GetPDGMass(); 
  }
  G4int       index = GetEnergyIndex(energy);
  G4double qeTotRat = GetNuMuQeTotRat(index, energy);

  G4ThreeVector dX = (lvX.vect()).unit();
  G4double eX   = lvX.e();  // excited nucleon
  G4double mX   = sqrt(massX2);
  G4double dP(0.), pX   = sqrt( eX*eX - mX*mX );
  G4double sumE = eX + rM;
  G4double a(0.), b(0.), c(0.), B(0.);

  if( qeTotRat > G4UniformRand() || mX <= fMt ) // || eX <= 1232.*MeV) // QE
  {  
    fString = false;

    if( fProton ) // pName == "nu_mu" ) 
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
    sumE = eX + rM; 
    G4double eTh = fMr+0.5*(fMr*fMr-mX*mX)/rM;

    if(eX <= eTh) // vmg, very rarely out of kinematics
    {
      theParticleChange.SetEnergyChange(energy);
      theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
      return &theParticleChange;
    } 
    B = sumE*sumE + rM*rM - fMr*fMr - pX*pX;
    a = 4.*(sumE*sumE - pX*pX);
    b = -4.*B*pX;
    c = 4.*sumE*sumE*rM*rM - B*B;
    G4double det2 = b*b-4.*a*c;
    if( det2 < 0.) det2 = 0.;
    dP = 0.5*(-b - sqrt(det2) )/a;
    pX -= dP;
    eX = sqrt( pX*pX + fMr*fMr );
    G4LorentzVector qeLV( pX*dX, eX );

    G4ParticleDefinition* qePart = G4ParticleTable::GetParticleTable()->
                                 FindParticle(fPDGencoding); 
 
    G4DynamicParticle* qeDyn = new G4DynamicParticle( qePart, qeLV);
    theParticleChange.AddSecondary(qeDyn); 
  
    G4double eRecoil = sqrt(rM*rM + dP*dP);
    G4ThreeVector vRecoil(dP*dX);
    G4LorentzVector lvTarg(vRecoil, eRecoil);

    if( eRecoil > 100.*MeV ) // add recoil nucleus
    {
      G4ParticleDefinition * recoilDef = 0;
      G4int Zr = recoil.GetZ_asInt();
      G4int Ar = recoil.GetA_asInt();

      if      ( Zr == 1 && Ar == 1 ) { recoilDef = G4Proton::Proton(); }
      else if ( Zr == 0 && Ar == 1 ) { recoilDef = G4Neutron::Neutron(); }
      else if ( Zr == 1 && Ar == 2 ) { recoilDef = G4Deuteron::Deuteron(); }
      else if ( Zr == 1 && Ar == 3 ) { recoilDef = G4Triton::Triton(); }
      else if ( Zr == 2 && Ar == 3 ) { recoilDef = G4He3::He3(); }
      else if ( Zr == 2 && Ar == 4 ) { recoilDef = G4Alpha::Alpha(); }
      else 
      {
        recoilDef = 
	G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon( Zr, Ar, 0.0 );
      }
      G4DynamicParticle * aSec = new G4DynamicParticle( recoilDef, lvTarg);
      theParticleChange.AddSecondary(aSec);
    } 
    else if( eRecoil > 0.0 ) 
    {
      theParticleChange.SetLocalEnergyDeposit( eRecoil );
    }
  }
  else if ( eX < 95000.*GeV ) // < 25.*GeV) //  < 95.*GeV ) // < 2.5*GeV ) //cluster decay
  {  
    if     (  fProton && pName == "nu_mu" )      qB =  1;
    else if(  fProton && pName == "anti_nu_mu" ) qB =  1;
    else if( !fProton && pName == "nu_mu" )      qB =  0;
    else if( !fProton && pName == "anti_nu_mu" ) qB =  0;

    // if( G4UniformRand() > 0.1 )
    {
      ClusterDecay( lvX, qB );
    }
    // else
    {
       pdgP = 111;

      if ( fQtransfer < 0.95*GeV )  // < 0.99*GeV )  //
      {
        // if( lvX.m() > mSum ) CoherentPion( lvX, pdgP, targetNucleus);
      }
    }
  }
  else // string
  {  
    return &theParticleChange;

  } 
  return &theParticleChange;
}


/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////
//
// sample x, then Q2

void G4NuMuNucleusNcModel::SampleLVkr(const G4HadProjectile & aTrack, G4Nucleus& targetNucleus)
{
  fBreak = false;
  G4int A = targetNucleus.GetA_asInt(), iTer(0), iTerMax(100); 
  G4int Z = targetNucleus.GetZ_asInt(); 
  G4double e3(0.), pMu2(0.), pX2(0.), nMom(0.), rM(0.), hM(0.), tM = targetNucleus.AtomicMass(A,Z);
  G4double cost(1.), sint(0.), phi(0.), muMom(0.); 
  G4ThreeVector eP, bst;
  const G4HadProjectile* aParticle = &aTrack;
  G4LorentzVector lvp1 = aParticle->Get4Momentum();
  nMom = NucleonMomentum( targetNucleus );

  if( A == 1 || nMom == 0. ) // hydrogen, no Fermi motion ???
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

      // if( e3 < sqrt(fW2) )  G4cout<<"energyX = "<<e3/GeV<<", fW = "<<sqrt(fW2)/GeV<<G4endl; // vmg ~10^-5 for NC
    
      pMu2 = fEmu*fEmu - fMnumu*fMnumu;
      pX2  = e3*e3 - fW2;

      fCosTheta  = fNuEnergy*fNuEnergy  + pMu2 - pX2;
      fCosTheta /= 2.*fNuEnergy*sqrt(pMu2);
      iTer++;
    }
    while( ( abs(fCosTheta) > 1. || fEmu < fMnumu ) && iTer < iTerMax );

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
    muMom = sqrt(fEmu*fEmu-fMnumu*fMnumu);
    eP *= muMom;
    fLVl = G4LorentzVector( eP, fEmu );

    fLVh = lvsum - fLVl;
    fLVt = G4LorentzVector( 0., 0., 0., 0. ); // no recoil
  }
  else // Fermi motion, Q2 in nucleon rest frame
  {
    G4ThreeVector nMomDir = nMom*G4RandomDirection();

    if( !f2p2h ) // 1p1h
    {
      G4Nucleus recoil(A-1,Z);
      rM = sqrt( recoil.AtomicMass(A-1,Z)*recoil.AtomicMass(A-1,Z) + nMom*nMom );
      hM = tM - rM;

      fLVt = G4LorentzVector( nMomDir, sqrt( rM*rM+nMom*nMom ) );
      fLVh = G4LorentzVector(-nMomDir, sqrt( hM*hM+nMom*nMom ) ); 
    }
    else // 2p2h
    {
      G4Nucleus recoil(A-2,Z-1);
      rM = recoil.AtomicMass(A-2,Z-1)+sqrt(nMom*nMom+fM1*fM1);
      hM = tM - rM;

      fLVt = G4LorentzVector( nMomDir, sqrt( rM*rM+nMom*nMom ) );
      fLVh = G4LorentzVector(-nMomDir, sqrt( hM*hM+nMom*nMom ) ); 
    }
    // G4cout<<hM<<", ";
    bst = fLVh.boostVector();

    lvp1.boost(-bst); // -> nucleon rest system, where Q2 transfer is ???

    fNuEnergy  = lvp1.e();
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

      // if(fEmu < 0.) G4cout<<"fEmu = "<<fEmu<<" hM = "<<hM<<G4endl;

      e3 = fNuEnergy + fM1 - fEmu;

      // if( e3 < sqrt(fW2) )  G4cout<<"energyX = "<<e3/GeV<<", fW = "<<sqrt(fW2)/GeV<<G4endl;
    
      pMu2 = fEmu*fEmu - fMnumu*fMnumu;
      pX2  = e3*e3 - fW2;

      fCosTheta  = fNuEnergy*fNuEnergy  + pMu2 - pX2;
      fCosTheta /= 2.*fNuEnergy*sqrt(pMu2);
      iTer++;
    }
    while( ( abs(fCosTheta) > 1. || fEmu < fMnumu ) && iTer < iTerMax );

    if( iTer >= iTerMax ) { fBreak = true; return; }

    if( abs(fCosTheta) > 1.) // vmg: due to big Q2/x values. To be improved ...
    { 
      G4cout<<"FM: fCosTheta = "<<fCosTheta<<", fEmu = "<<fEmu<<G4endl;
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
    muMom = sqrt(fEmu*fEmu-fMnumu*fMnumu);
    eP *= muMom;
    fLVl = G4LorentzVector( eP, fEmu );
    fLVh = lvsum - fLVl;
    // back to lab system
    fLVl.boost(bst);
    fLVh.boost(bst);
  }
  //G4cout<<iTer<<", "<<fBreak<<"; ";
}

//////////////////////////////////////

G4double G4NuMuNucleusNcModel::SampleXkr(G4double energy)
{
  G4int i(0), nBin(50);
  G4double xx(0.), prob = G4UniformRand();

  for( i = 0; i < nBin; ++i ) 
  {
    if( energy <= fNuMuEnergyLogVector[i] ) break;
  }
  if( i <= 0)          // E-edge
  {
    fEindex = 0;
    xx = GetXkr( 0, prob);  
  }
  else if ( i >= nBin-1) 
  {
    fEindex = nBin-1;  
    xx = GetXkr( nBin-1, prob); 
  }
  else
  {
    fEindex = i;
    G4double x1 = GetXkr(i-1,prob);
    G4double x2 = GetXkr(i,prob);

    G4double e1 = G4Log(fNuMuEnergyLogVector[i-1]);
    G4double e2 = G4Log(fNuMuEnergyLogVector[i]);
    G4double e  = G4Log(energy);

    if( e2 <= e1) xx = x1 + G4UniformRand()*(x2-x1);
    else          xx = x1 + (e-e1)*(x2-x1)/(e2-e1);  // lin in energy log-scale
  }
  return xx;
}

//////////////////////////////////////////////
//
// sample X according to prob (xmin,1) at a given energy index iEnergy

G4double G4NuMuNucleusNcModel::GetXkr(G4int iEnergy, G4double prob)
{
  G4int i(0), nBin=50; 
  G4double xx(0.);

  for( i = 0; i < nBin; ++i ) 
  {
    if( prob <= fNuMuXdistrKR[iEnergy][i] ) 
      break;
  }
  if(i <= 0 )  // X-edge
  {
    fXindex = 0;
    xx = fNuMuXarrayKR[iEnergy][0];
  }
  if ( i >= nBin ) 
  {
    fXindex = nBin;
    xx = fNuMuXarrayKR[iEnergy][nBin];
  }  
  else
  {
    fXindex = i;
    G4double x1 = fNuMuXarrayKR[iEnergy][i];
    G4double x2 = fNuMuXarrayKR[iEnergy][i+1];

    G4double p1 = 0.;

    if( i > 0 ) p1 = fNuMuXdistrKR[iEnergy][i-1];

    G4double p2 = fNuMuXdistrKR[iEnergy][i];  

    if( p2 <= p1 ) xx = x1 + G4UniformRand()*(x2-x1);
    else           xx = x1 + (prob-p1)*(x2-x1)/(p2-p1);
  }
  return xx;
}

//////////////////////////////////////
//
// Sample fQtransfer at a given Enu and fX

G4double G4NuMuNucleusNcModel::SampleQkr( G4double energy, G4double xx)
{
  G4int nBin(50), iE=fEindex, jX=fXindex;
  G4double qq(0.), qq1(0.), qq2(0.);
  G4double prob = G4UniformRand();

  // first E

  if( iE <= 0 )          
  {
    qq1 = GetQkr( 0, jX, prob);  
  }
  else if ( iE >= nBin-1) 
  {
    qq1 = GetQkr( nBin-1, jX, prob); 
  }
  else
  {
    G4double q1 = GetQkr(iE-1,jX, prob);
    G4double q2 = GetQkr(iE,jX, prob);

    G4double e1 = G4Log(fNuMuEnergyLogVector[iE-1]);
    G4double e2 = G4Log(fNuMuEnergyLogVector[iE]);
    G4double e  = G4Log(energy);

    if( e2 <= e1) qq1 = q1 + G4UniformRand()*(q2-q1);
    else          qq1 = q1 + (e-e1)*(q2-q1)/(e2-e1);  // lin in energy log-scale
  }

  // then X

  if( jX <= 0 )          
  {
    qq2 = GetQkr( iE, 0, prob);  
  }
  else if ( jX >= nBin) 
  {
    qq2 = GetQkr( iE, nBin, prob); 
  }
  else
  {
    G4double q1 = GetQkr(iE,jX-1, prob);
    G4double q2 = GetQkr(iE,jX, prob);

    G4double e1 = G4Log(fNuMuXarrayKR[iE][jX-1]);
    G4double e2 = G4Log(fNuMuXarrayKR[iE][jX]);
    G4double e  = G4Log(xx);

    if( e2 <= e1) qq2 = q1 + G4UniformRand()*(q2-q1);
    else          qq2 = q1 + (e-e1)*(q2-q1)/(e2-e1);  // lin in energy log-scale
  }
  qq = 0.5*(qq1+qq2);

  return qq;
}

//////////////////////////////////////////////
//
// sample Q according to prob (qmin,qmax) at a given energy index iE and X index jX

G4double G4NuMuNucleusNcModel::GetQkr( G4int iE, G4int jX, G4double prob )
{
  G4int i(0), nBin=50; 
  G4double qq(0.);

  for( i = 0; i < nBin; ++i ) 
  {
    if( prob <= fNuMuQdistrKR[iE][jX][i] ) 
      break;
  }
  if(i <= 0 )  // Q-edge
  {
    fQindex = 0;
    qq = fNuMuQarrayKR[iE][jX][0];
  }
  if ( i >= nBin ) 
  {
    fQindex = nBin;
    qq = fNuMuQarrayKR[iE][jX][nBin];
  }  
  else
  {
    fQindex = i;
    G4double q1 = fNuMuQarrayKR[iE][jX][i];
    G4double q2 = fNuMuQarrayKR[iE][jX][i+1];

    G4double p1 = 0.;

    if( i > 0 ) p1 = fNuMuQdistrKR[iE][jX][i-1];

    G4double p2 = fNuMuQdistrKR[iE][jX][i];  

    if( p2 <= p1 ) qq = q1 + G4UniformRand()*(q2-q1);
    else           qq = q1 + (prob-p1)*(q2-q1)/(p2-p1);
  }
  return qq;
}

//
//
///////////////////////////
