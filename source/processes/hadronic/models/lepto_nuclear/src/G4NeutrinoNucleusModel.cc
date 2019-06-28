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
// $Id: G4NeutrinoNucleusModel.cc 91806 2015-08-06 12:20:45Z gcosmo $
//
// Geant4 Header : G4NeutrinoNucleusModel
//
// Author : V.Grichine 12.2.19
//  

#include "G4NeutrinoNucleusModel.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"

//#include "G4Integrator.hh"
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"

#include "G4CascadeInterface.hh"
// #include "G4BinaryCascade.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4FTFModel.hh"
// #include "G4BinaryCascade.hh"
#include "G4HadFinalState.hh"
#include "G4HadSecondary.hh"
#include "G4HadronicInteractionRegistry.hh"
// #include "G4INCLXXInterface.hh"
#include "G4KineticTrack.hh"
#include "G4DecayKineticTracks.hh"
#include "G4KineticTrackVector.hh"
#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"

// #include "G4QGSModel.hh"
// #include "G4QGSMFragmentation.hh"
// #include "G4QGSParticipants.hh"


#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Nucleus.hh"
#include "G4LorentzVector.hh"

using namespace std;
using namespace CLHEP;

const G4int G4NeutrinoNucleusModel::fResNumber = 6;

const G4double G4NeutrinoNucleusModel::fResMass[6] = // [fResNumber] = 
  {2190., 1920., 1700., 1600., 1440., 1232. };

const G4int G4NeutrinoNucleusModel::fClustNumber = 4;

const G4double G4NeutrinoNucleusModel::fMesMass[4] = {1260., 980., 770., 139.57};
const G4int    G4NeutrinoNucleusModel::fMesPDG[4]  = {20213, 9000211, 213, 211};

// const G4double G4NeutrinoNucleusModel::fBarMass[4] = {1905., 1600., 1232., 939.57};
// const G4int    G4NeutrinoNucleusModel::fBarPDG[4]  = {2226, 32224, 2224, 2212};

const G4double G4NeutrinoNucleusModel::fBarMass[4] = {1700., 1600., 1232., 939.57};
const G4int    G4NeutrinoNucleusModel::fBarPDG[4]  = {12224, 32224, 2224, 2212};


G4NeutrinoNucleusModel::G4NeutrinoNucleusModel(const G4String& name) 
  : G4HadronicInteraction(name)
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( 100.*TeV );
  SetMinEnergy(1.e-6*eV); 

  fNbin = 50;
  fEindex = fXindex = 0;
  fOnePionIndex = 58;
  fIndex = 50;
  fCascade = fString = fProton = f2p2h = false;

  fNuEnergy  = fQ2  = fQtransfer = fXsample = fDp = 0.;
  fCosTheta = fCosThetaPi = 1.;
  fEmuPi = fW2 = fW2pi = 0.;

  fMu = 105.6583745*MeV;
  fMpi = 139.57018*MeV;
  fM1 = 939.5654133*MeV;  // for nu_mu -> mu-,  and n -> p
  fM2 = 938.2720813*MeV;

  fEmu = fMu;
  fEx = fM1;
  fMr = 1232.*MeV;
  fMt = fM2; // threshold for N*-diffraction

  fMinNuEnergy = GetMinNuMuEnergy();

  fLVh = G4LorentzVector(0.,0.,0.,0.);
  fLVl = G4LorentzVector(0.,0.,0.,0.);
  fLVt = G4LorentzVector(0.,0.,0.,0.);
  fLVcpi = G4LorentzVector(0.,0.,0.,0.);

  // theMuonMinus = G4MuonMinus::MuonMinus();
  // theMuonPlus = G4MuonPlus::MuonPlus();

  // PDG2016: sin^2 theta Weinberg

  fSin2tW = 0.23129; // 0.2312;

  fCutEnergy = 0.; // default value

  fPDGencoding = 0; // unphysical as default

// reuse existing pre-compound model
/*
  G4GeneratorPrecompoundInterface* precoInterface = new G4GeneratorPrecompoundInterface();

  G4HadronicInteraction* p = G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");

  fPrecoModel = static_cast<G4VPreCompoundModel*>(p);

  if(!fPrecoModel) fPrecoModel = new G4PreCompoundModel(); 

  precoInterface->SetDeExcitation(fPrecoModel);

  // binary with fPrecoModel

  // theBinary = new G4BinaryCascade(fPrecoModel);

  // INCLXX with fPrecoModel

  // theINCLXX = new G4INCLXXInterface(fPrecoModel);

  // Build Bertini model

  theBertini = new G4CascadeInterface();

  // FTFP string model

  theFTFP = new G4TheoFSGenerator();
  theFTFP->SetTransport(precoInterface);
  theFragmentation = new G4LundStringFragmentation();
  theStringDecay = new G4ExcitedStringDecay(theFragmentation);
  G4FTFModel* theStringModel = new G4FTFModel();
  theStringModel->SetFragmentationModel(theStringDecay);
  theFTFP->SetHighEnergyGenerator(theStringModel);

  //  QGSP string model


  theQGSP = new G4TheoFSGenerator("QGSP");

  G4QGSModel< G4QGSParticipants >* stringModel = new G4QGSModel< G4QGSParticipants >;
  G4ExcitedStringDecay* stringDecay = new G4ExcitedStringDecay(new G4QGSMFragmentation);
  stringModel->SetFragmentationModel(stringDecay);

  // theCascade = new G4GeneratorPrecompoundInterface();

  theQGSP->SetTransport(precoInterface);
  theQGSP->SetHighEnergyGenerator(stringModel);
*/
  fRecoil = nullptr;
     
}


G4NeutrinoNucleusModel::~G4NeutrinoNucleusModel()
{}


void G4NeutrinoNucleusModel::ModelDescription(std::ostream& outFile) const
{

    outFile << "G4NeutrinoNucleusModel is a neutrino-nucleus general\n"
            << "model which uses the standard model \n"
            << "transfer parameterization.  The model is fully relativistic\n";

}

/////////////////////////////////////////////////////////

G4bool G4NeutrinoNucleusModel::IsApplicable(const G4HadProjectile & aPart, 
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
/*
/////////////////////////////////////////// ClusterDecay ////////////////////////////////////////////////////////////
//
//

G4HadFinalState* G4NeutrinoNucleusModel::ApplyYourself(
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

    // muMom = sqrt(fEmuPi*fEmuPi-fMu*fMu);
    muMom = sqrt(fEmu*fEmu-fMu*fMu);

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

    if(  pName == "nu_mu" )         aLept = new G4DynamicParticle( theMuonMinus, lv2 );  
    else if( pName == "anti_nu_mu") aLept = new G4DynamicParticle( theMuonPlus,  lv2 );

    if( pName == "nu_mu" ) pdgP =  211;
    else                   pdgP = -211;

    G4double eCut = fMpi + 0.5*(fMpi*fMpi-massX2)/mTarg; // massX -> fMpi

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

    muMom = sqrt(fEmu*fEmu-fMu*fMu);

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

    if(  pName == "nu_mu" )         aLept = new G4DynamicParticle( theMuonMinus, lv2 );  
    else if( pName == "anti_nu_mu") aLept = new G4DynamicParticle( theMuonPlus,  lv2 );
    
    theParticleChange.AddSecondary( aLept );
  }

  // hadron part

  fRecoil  = nullptr;
  fCascade = false;
  fString  = false;
  
  if( A == 1 )
  {
    if( pName == "nu_mu" ) qB = 2;
    else                   qB = 0;

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

    if( pName == "nu_mu" ) // (++) state -> p + pi+
    { 
      fMt = G4ParticleTable::GetParticleTable()->FindParticle(2212)->GetPDGMass()
          + G4ParticleTable::GetParticleTable()->FindParticle(211)->GetPDGMass();
    }
    else // (0) state -> p + pi-, n + pi0
    {
      fMt = G4ParticleTable::GetParticleTable()->FindParticle(2212)->GetPDGMass()
          + G4ParticleTable::GetParticleTable()->FindParticle(-211)->GetPDGMass();
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
      fMt = G4ParticleTable::GetParticleTable()->FindParticle(2112)->GetPDGMass()
          + G4ParticleTable::GetParticleTable()->FindParticle(-211)->GetPDGMass();
    } 
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

    if( pName == "nu_mu" ) 
    {  
      fPDGencoding = 2212;
      fMr =  proton_mass_c2;
      recoil = G4Nucleus(A-1,Z);
      fRecoil = &recoil;
      rM = recoil.AtomicMass(A-1,Z);
    } 
    else // if( pName == "anti_nu_mu" ) 
    {  
      fPDGencoding = 2112;
      fMr =   G4ParticleTable::GetParticleTable()->
	FindParticle(fPDGencoding)->GetPDGMass(); // 939.5654133*MeV;
      recoil = G4Nucleus(A-1,Z-1);
      fRecoil = &recoil;
      rM = recoil.AtomicMass(A-1,Z-1);
    } 
    sumE = eX + rM;   
    B = sumE*sumE + rM*rM - fMr*fMr - pX*pX;
    a = 4.*(sumE*sumE - pX*pX);
    b = -4.*B*pX;
    c = 4.*sumE*sumE*rM*rM - B*B;
    dP = 0.5*(-b - sqrt(b*b-4.*a*c) )/a;
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
  else if ( eX < 25.*GeV) //  < 95000.*GeV ) // < 95.*GeV ) // < 2.5*GeV ) //cluster decay
  {  
    if     (  fProton && pName == "nu_mu" )      qB =  2;
    else if(  fProton && pName == "anti_nu_mu" ) qB =  0;
    else if( !fProton && pName == "nu_mu" )      qB =  1;
    else if( !fProton && pName == "anti_nu_mu" ) qB = -1;

    // if( G4UniformRand() > 0.1 )
    {
      ClusterDecay( lvX, qB );
    }
    // else
    {
      if( pName == "nu_mu" ) pdgP =  211;
      else                   pdgP = -211;

      if ( fQtransfer < 0.95*GeV )  // < 0.99*GeV )  //
      {
        // if( lvX.m() > mSum ) CoherentPion( lvX, pdgP, targetNucleus);
      }
    }
  }
  else // string
  {  
    return &theParticleChange;

    fString = true;
     
    if( pName == "nu_mu" ) 
    {  
      fPDGencoding = 2212;
      fMr =  proton_mass_c2;
      recoil = G4Nucleus(A-1,Z);
      fRecoil = &recoil;
    }
    else // if( pName == "anti_nu_mu" ) 
    {  
      fPDGencoding = 2112;
      fMr =  939.5654133*MeV;
      recoil = G4Nucleus(A-1,Z-1);
      fRecoil = &recoil;
    }
    pX = sqrt( eX*eX - fMr*fMr );
    G4LorentzVector qeLV( pX*dX, eX );

    G4ParticleDefinition* qePart = G4ParticleTable::GetParticleTable()->
                                 FindParticle(fPDGencoding); 
 
    G4DynamicParticle qeDyn( qePart, qeLV);
    G4HadProjectile projectile(qeDyn);

    // G4HadFinalState* hfs = theFTFP->ApplyYourself(projectile, recoil);
    G4HadFinalState* hfs = theQGSP->ApplyYourself(projectile, recoil);

    theParticleChange.AddSecondaries( hfs );
  } 
  return &theParticleChange;
}


/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////
//
// sample x, then Q2

void G4NeutrinoNucleusModel::SampleLVkr(const G4HadProjectile & aTrack, G4Nucleus& targetNucleus)
{
  fBreak = false;
  G4int A = targetNucleus.GetA_asInt(), iTer(0), iTerMax(20); 
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

      if( e3 < sqrt(fW2) )  G4cout<<"energyX = "<<e3/GeV<<", fW = "<<sqrt(fW2)/GeV<<G4endl;
    
      pMu2 = fEmu*fEmu - fMu*fMu;
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
    G4ThreeVector nMomDir = nMom*G4RandomDirection();

    if( !f2p2h ) // 1p1h
    {
      G4Nucleus recoil(A-1,Z);
      rM = sqrt( recoil.AtomicMass(A-1,Z)*recoil.AtomicMass(A-1,Z) + nMom*nMom );
      hM = tM - rM;

      fLVt = G4LorentzVector( nMomDir, rM );
      fLVh = G4LorentzVector(-nMomDir, hM); 
    }
    else // 2p2h
    {
      G4Nucleus recoil(A-2,Z-1);
      rM = recoil.AtomicMass(A-2,Z-1)+sqrt(nMom*nMom+fM1*fM1);
      hM = tM - rM;

      fLVt = G4LorentzVector( nMomDir, rM );
      fLVh = G4LorentzVector(-nMomDir, hM); 
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
    
      pMu2 = fEmu*fEmu - fMu*fMu;
      pX2  = e3*e3 - fW2;

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
    // back to lab system
    fLVl.boost(bst);
    fLVh.boost(bst);
  }
  //G4cout<<iTer<<", "<<fBreak<<"; ";
}

//////////////////////////////////////

G4double G4NeutrinoNucleusModel::SampleXkr(G4double energy)
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
  else if ( i >= nBin) 
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

G4double G4NeutrinoNucleusModel::GetXkr(G4int iEnergy, G4double prob)
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

G4double G4NeutrinoNucleusModel::SampleQkr( G4double energy, G4double xx)
{
  G4int nBin(50), iE=fEindex, jX=fXindex;
  G4double qq(0.), qq1(0.), qq2(0.);
  G4double prob = G4UniformRand();

  // first E

  if( iE <= 0 )          
  {
    qq1 = GetQkr( 0, jX, prob);  
  }
  else if ( iE >= nBin) 
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
  else if ( iE >= nBin) 
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

G4double G4NeutrinoNucleusModel::GetQkr( G4int iE, G4int jX, G4double prob )
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
    fXindex = 0;
    qq = fNuMuQarrayKR[iE][jX][0];
  }
  if ( i >= nBin ) 
  {
    fXindex = nBin;
    qq = fNuMuQarrayKR[iE][jX][nBin];
  }  
  else
  {
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

*/

///////////////////////////////////////////////////////////
//
// Final meson to theParticleChange

void G4NeutrinoNucleusModel::FinalMeson( G4LorentzVector & lvM, G4int, G4int pdgM) // qM
{
  G4int pdg = pdgM;
  // if      ( qM ==  0 ) pdg = pdgM - 100;
  // else if ( qM == -1 ) pdg = -pdgM;

  if( pdg == 211 || pdg == -211 || pdg == 111) // pions
  {
    G4ParticleDefinition* pd2 = G4ParticleTable::GetParticleTable()->FindParticle(pdg); 
    G4DynamicParticle*    dp2 = new G4DynamicParticle( pd2, lvM);
    theParticleChange.AddSecondary( dp2 );
  }
  else // meson resonances
  {
    G4ParticleDefinition* rePart = G4ParticleTable::GetParticleTable()->
      FindParticle(pdg); 
    G4KineticTrack ddkt( rePart, 0., G4ThreeVector(0.,0.,0.), lvM);
    G4KineticTrackVector* ddktv = ddkt.Decay();

    G4DecayKineticTracks decay( ddktv );

    for( unsigned int i = 0; i < ddktv->size(); i++ ) // add products to partchange
    {
      G4DynamicParticle * aNew = 
      new G4DynamicParticle( ddktv->operator[](i)->GetDefinition(),
                             ddktv->operator[](i)->Get4Momentum());

      // G4cout<<"       "<<i<<", "<<aNew->GetDefinition()->GetParticleName()<<", "<<aNew->Get4Momentum()<<G4endl;

      theParticleChange.AddSecondary( aNew );
      delete ddktv->operator[](i);
    }
    delete ddktv;
  }
}

////////////////////////////////////////////////////////
//
// Final barion to theParticleChange, and recoil nucleus treatment

void G4NeutrinoNucleusModel::FinalBarion( G4LorentzVector & lvB, G4int, G4int pdgB) // qB
{
  G4int A(0), Z(0), pdg = pdgB;

  // if      ( qB ==  1 ) pdg = pdgB - 10;
  // else if ( qB ==  0 ) pdg = pdgB - 110;
  // else if ( qB == -1 ) pdg = pdgB - 1110;

  if( pdg == 2212 || pdg == 2112) fMr = G4ParticleTable::GetParticleTable()->FindParticle(pdg)->GetPDGMass();
  else fMr = lvB.m();
  G4double eX = lvB.e();
  G4double rM(0.), mX = lvB.m();
  G4ThreeVector dX = (lvB.vect()).unit();
  G4double pX = sqrt(eX*eX-mX*mX);

  if( fRecoil )
  {
    Z  = fRecoil->GetZ_asInt();
    A  = fRecoil->GetA_asInt();
    rM = fRecoil->AtomicMass(A,Z); //->AtomicMass(); // 
  }
  else // A=0 nu+p
  {
    A = 0;
    Z = 1;
    rM = electron_mass_c2;
  }
  // G4cout<<A<<", ";

  G4double sumE = eX + rM;   
  G4double B = sumE*sumE + rM*rM - fMr*fMr - pX*pX;
  G4double a = 4.*(sumE*sumE - pX*pX);
  G4double b = -4.*B*pX;
  G4double c = 4.*sumE*sumE*rM*rM - B*B;
  G4double dP = 0.5*(-b - sqrt(b*b-4.*a*c) )/a;

  fDp = dP;

  pX -= dP;

  // if( A == 0 ) G4cout<<pX/MeV<<", ";

  eX = sqrt( pX*pX + fMr*fMr );
  G4LorentzVector lvN( pX*dX, eX );

  if( pdg == 2212 || pdg == 2112) // nucleons mX >= fMr, dP >= 0
  {
    G4ParticleDefinition* pd2 = G4ParticleTable::GetParticleTable()->FindParticle(pdg); 
    G4DynamicParticle*    dp2 = new G4DynamicParticle( pd2, lvN);
    theParticleChange.AddSecondary( dp2 );

  }
  else // delta resonances
  {
    G4ParticleDefinition* rePart = G4ParticleTable::GetParticleTable()->FindParticle(pdg); 
    G4KineticTrack ddkt( rePart, 0., G4ThreeVector(0.,0.,0.), lvN);
    G4KineticTrackVector* ddktv = ddkt.Decay();

    G4DecayKineticTracks decay( ddktv );

    for( unsigned int i = 0; i < ddktv->size(); i++ ) // add products to partchange
    {
      G4DynamicParticle * aNew = 
      new G4DynamicParticle( ddktv->operator[](i)->GetDefinition(),
                             ddktv->operator[](i)->Get4Momentum());

      // G4cout<<"       "<<i<<", "<<aNew->GetDefinition()->GetParticleName()<<", "<<aNew->Get4Momentum()<<G4endl;

      theParticleChange.AddSecondary( aNew );
      delete ddktv->operator[](i);
    }
    delete ddktv;
  }
  // recoil nucleus

  G4double eRecoil = sqrt( rM*rM + dP*dP );
  fTr = eRecoil - rM;
  G4ThreeVector vRecoil(dP*dX);
  G4LorentzVector lvTarg(vRecoil, eRecoil);

  G4ParticleDefinition* recoilDef = 0;

  // if( G4UniformRand() > 0.5 )
  if( G4UniformRand() >= 0.0 )
  {
    if( fTr > 100.*MeV && A > 0 ) // add recoil nucleus
    {

      if      ( Z == 1 && A == 1 ) { recoilDef = G4Proton::Proton(); }
      else if ( Z == 0 && A == 1 ) { recoilDef = G4Neutron::Neutron(); }
      else if ( Z == 1 && A == 0 ) { recoilDef = G4Positron::Positron(); } // dP to positron, if nu+p
      else if ( Z == 1 && A == 2 ) { recoilDef = G4Deuteron::Deuteron(); }
      else if ( Z == 1 && A == 3 ) { recoilDef = G4Triton::Triton(); }
      else if ( Z == 2 && A == 3 ) { recoilDef = G4He3::He3(); }
      else if ( Z == 2 && A == 4 ) { recoilDef = G4Alpha::Alpha(); }
      else 
      {
        recoilDef = 
	G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon( Z, A, 0.0 );
      }
      G4DynamicParticle * aSec = new G4DynamicParticle( recoilDef, lvTarg);
      theParticleChange.AddSecondary(aSec);
    } 
    else if( eRecoil > 0.0 ) 
    {
      if ( A > 0 ) theParticleChange.SetLocalEnergyDeposit( eRecoil );
      else theParticleChange.SetLocalEnergyDeposit( dP ); // recoil momentum as energy deposition
    }
  }
  else if( A > 0)
  {
    G4ThreeVector bst(0.,0.,0.);
    G4LorentzVector lvR( bst, eRecoil);
    G4Fragment* fragment = new G4Fragment(A,Z,lvR);
    fragment->SetNumberOfHoles(1);

    /*
    G4VPreCompaundModel* dexcite = fPre;

    G4ReactionProductVector* products = fPrecoModel->DeExcite(*fragment); 

    G4ReactionProductVector::iterator iter;

    for(iter = products->begin(); iter != products->end(); ++iter)
    {
      G4DynamicParticle * aNewDP =
                    new G4DynamicParticle((*iter)->GetDefinition(),
                            (*iter)->GetTotalEnergy(),
                            (*iter)->GetMomentum());
      G4HadSecondary aNew = G4HadSecondary(aNewDP);

      G4double time=(*iter)->GetFormationTime();

      if(time < 0.0) { time = 0.0; }

      aNew.SetTime(time);// (timePrimary + time);
      aNew.SetCreatorModelType((*iter)->GetCreatorModel());

      theParticleChange.AddSecondary(aNew);
    }
    */
    delete fragment;
    fragment = nullptr; 
  }
  else // 
  {
    theParticleChange.SetLocalEnergyDeposit( eRecoil );
    /*
    recoilDef = G4Positron::Positron();
    G4DynamicParticle * aSec = new G4DynamicParticle( recoilDef, lvTarg);
    theParticleChange.AddSecondary(aSec);
    */
  }
}

///////////////////////////////////////////
//
// Fragmentation of lvX directly to pion and recoil nucleus (A,Z)

void G4NeutrinoNucleusModel::CoherentPion( G4LorentzVector & lvP, G4int pdgP, G4Nucleus & targetNucleus)
{
  G4int A(0), Z(0), pdg = pdgP;
  fLVcpi = G4LorentzVector(0.,0.,0.,0.);

  G4double  rM(0.), mN(938.), mI(0.); 

  mN = G4ParticleTable::GetParticleTable()->FindParticle(2212)->GetPDGMass(); // *0.85; // *0.9; // 

  // mN = 1.*139.57 + G4UniformRand()*(938. - 1.*139.57);

  G4ThreeVector vN = lvP.boostVector(), bst(0.,0.,0.);
  //  G4double gN = lvP.e()/lvP.m();
  // G4LorentzVector  lvNu(vN*gN*mN, mN*gN);
  G4LorentzVector  lvNu(bst, mN);

  // lvP = lvP - lvNu; // already 1pi

  // G4cout<<vN-lvP.boostVector()<<", ";

  Z  = targetNucleus.GetZ_asInt();
  A  = targetNucleus.GetA_asInt();
  rM = targetNucleus.AtomicMass(A,Z); //->AtomicMass(); // 

  // G4cout<<rM<<", ";
  // G4cout<<A<<", ";
  
  if( A == 1 ) 
  {
    // bst = lvNu.boostVector();
    mI = 0.;
  }
  else
  {
    G4Nucleus targ(A-1,Z);
    mI = targ.AtomicMass(A-1,Z);
    G4LorentzVector lvTar(bst,rM); 
    lvNu = lvNu + lvTar;
    // bst = lvNu.boostVector();  
    bst = fLVt.boostVector();  
    lvP.boost(-bst);
  }
  fMr = G4ParticleTable::GetParticleTable()->FindParticle(pdg)->GetPDGMass();
  G4double eX = lvP.e();
  G4double mX = lvP.m();
  // G4cout<<mX-fMr<<", ";
  G4ThreeVector dX = (lvP.vect()).unit();
  // G4cout<<dX<<", ";
  G4double pX = sqrt(eX*eX-mX*mX);
  // G4cout<<pX<<", ";
  G4double sumE = eX + rM;   
  G4double B = sumE*sumE + rM*rM - fMr*fMr - pX*pX;
  G4double a = 4.*(sumE*sumE - pX*pX);
  G4double b = -4.*B*pX;
  G4double c = 4.*sumE*sumE*rM*rM - B*B;
  G4double dP = 0.5*(-b - sqrt(b*b-4.*a*c) )/a;

  dP = FinalMomentum( mI, rM, fMr, lvP);

  // G4cout<<dP<<", ";
  pX -= dP;
  eX = sqrt( pX*pX + fMr*fMr );
  G4LorentzVector lvN( pX*dX, eX );

  fLVcpi = lvN;

  if( A > 1 ) lvN.boost(bst);

  G4ParticleDefinition* pd2 = G4ParticleTable::GetParticleTable()->FindParticle(pdg); 
  G4DynamicParticle*    dp2 = new G4DynamicParticle( pd2, lvN);
  theParticleChange.AddSecondary( dp2 );

  // recoil nucleus

  G4double eRecoil = sqrt( rM*rM + dP*dP );
  G4ThreeVector vRecoil(dP*dX);
  G4LorentzVector lvTarg(vRecoil, eRecoil);
  // lvTarg.boost(bst);

  // G4LorentzVector lvSum = lvN+lvTarg; G4cout<<lvSum.m()/GeV<<", ";

  if( eRecoil > 0.*MeV ) //100.*MeV ) // add recoil nucleus
  {
    G4ParticleDefinition * recoilDef = 0;

    if      ( Z == 1 && A == 1 ) { recoilDef = G4Proton::Proton(); }
    else if ( Z == 0 && A == 1 ) { recoilDef = G4Neutron::Neutron(); }
    else if ( Z == 1 && A == 0 ) { recoilDef = G4Positron::Positron(); } // dP to positron, if nu+p
    else if ( Z == 1 && A == 2 ) { recoilDef = G4Deuteron::Deuteron(); }
    else if ( Z == 1 && A == 3 ) { recoilDef = G4Triton::Triton(); }
    else if ( Z == 2 && A == 3 ) { recoilDef = G4He3::He3(); }
    else if ( Z == 2 && A == 4 ) { recoilDef = G4Alpha::Alpha(); }
    else 
    {
        recoilDef = 
	G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon( Z, A, 0.0 );
    }
    G4DynamicParticle * aSec = new G4DynamicParticle( recoilDef, lvTarg);
    theParticleChange.AddSecondary(aSec);
  } 
  else if( eRecoil > 0.0 ) 
  {
      theParticleChange.SetLocalEnergyDeposit( eRecoil );
  }
}

////////////////////////////////////////////////////////////
//
// Excited barion decay to meson and barion, 
// mass distributions and charge exchange are free parameters 

void G4NeutrinoNucleusModel::ClusterDecay( G4LorentzVector & lvX, G4int qX)
{
  G4bool   finB = false;
  G4int    pdgB(0), i(0), qM(0), qB(0); //   pdgM(0),
  G4double mM(0.), mB(0.), eM(0.), eB(0.), pM(0.), pB(0.);
  G4double mm1(0.), mm22(0.), M1(0.), M2(0.), mX(0.);

  mX = lvX.m();

  G4double mN  = G4ParticleTable::GetParticleTable()->FindParticle(2112)->GetPDGMass();
  G4double mPi = G4ParticleTable::GetParticleTable()->FindParticle(211)->GetPDGMass();

  //  G4double deltaM =  1.*MeV; // 30.*MeV; //  10.*MeV; //  100.*MeV; // 20.*MeV; // 
  G4double deltaMr[4] =  { 0.*MeV, 0.*MeV, 100.*MeV, 0.*MeV}; 

  G4ThreeVector   dir(0.,0.,0.);
  G4ThreeVector   bst(0.,0.,0.);
  G4LorentzVector lvM(0.,0.,0.,0.);
  G4LorentzVector lvB(0.,0.,0.,0.);

  for( i = 0; i < fClustNumber; ++i) // check resonance
  {
    if( mX >= fBarMass[i] )
    {
        pdgB = fBarPDG[i];
        // mB = G4ParticleTable::GetParticleTable()->FindParticle(pdgB)->GetPDGMass();
        break;
    }
  }
  if( i == fClustNumber || i == fClustNumber-1 ) // low mass, p || n
  {
    if     ( qX == 2 || qX ==  0) { pdgB = 2212; qB = 1;} // p for 2, 0
    
    else if( qX == 1 || qX == -1) { pdgB = 2112; qB = 0;} // n for 1, -1

    return FinalBarion( lvX, qB, pdgB);     
  }
  else if( mX < fBarMass[i] + deltaMr[i]  || mX < mN + mPi ) 
  {
    finB = true; // final barion -> out

    if     ( qX == 1  && pdgB != 2212)  pdgB = pdgB - 10;
    else if( qX == 0  && pdgB != 2212)  pdgB = pdgB - 110;
    else if( qX == 0  && pdgB == 2212)  pdgB = pdgB - 100;

    if( finB ) return FinalBarion( lvX, qX, pdgB ); // out
  }
  // no barion resonance, try 1->2 decay in COM frame

  // try meson mass

  mm1 = mPi + 1.*MeV;     // pi+
  mm22 = mX - mN; // mX-n

  if( mm22 <= mm1 ) // out with p or n
  {
    if( qX == 2 || qX == 0)       { pdgB = 2212; qB = 1;} // p
    else if( qX == 1 || qX == -1) { pdgB = 2112; qB = 0;} // n

    return FinalBarion(lvX, qB, pdgB);
  }
  else // try decay -> meson(cluster) + barion(cluster)
  {
    // G4double sigmaM = 50.*MeV; // 100.*MeV; // 200.*MeV; // 400.*MeV; // 800.*MeV; //
    G4double rand   = G4UniformRand();

    // mM = mm1*mm22/( mm1 + rand*(mm22 - mm1) );
    // mM = mm1*mm22/sqrt( mm1*mm1 + rand*(mm22*mm22 - mm1*mm1) );
    // mM = -sigmaM*log( (1.- rand)*exp(-mm22/sigmaM) + rand*exp(-mm1/sigmaM) );
    mM = mm1 + rand*(mm22-mm1);


    for( i = 0; i < fClustNumber; ++i)
    {
      if( mM >= fMesMass[i] )
      {
        // pdgM = fMesPDG[i];
        // mM = G4ParticleTable::GetParticleTable()->FindParticle(pdgM)->GetPDGMass();
        break;
      }
    }
    M1 = G4ParticleTable::GetParticleTable()->FindParticle(2112)->GetPDGMass()+2.*MeV; // n
    M2 = mX - mM;

    if( M2 <= M1 ) // 
    {
      if     ( qX == 2 || qX ==  0) { pdgB = 2212; qB = 1;} // p
      else if( qX == 1 || qX == -1) { pdgB = 2112; qB = 0;} // n

      return FinalBarion(lvX, qB, pdgB);     
    }
    mB = M1 + G4UniformRand()*(M2-M1);
    // mB = -sigmaM*log( (1.- rand)*exp(-M2/sigmaM) + rand*exp(-M1/sigmaM) );


    dir = G4RandomDirection(); // ???
    bst = lvX.boostVector();

    eM  = 0.5*(mX*mX + mM*mM - mB*mB)/mX;
    pM  = sqrt(eM*eM - mM*mM);
    lvM = G4LorentzVector( pM*dir, eM);
    lvM.boost(bst);

    eB  = 0.5*(mX*mX + mB*mB - mM*mM)/mX;
    pB  = sqrt(eB*eB - mB*mB);
    lvB = G4LorentzVector(-pB*dir, eB);
    lvB.boost(bst);

    // G4cout<<mM<<"/"<<mB<<", ";

    // charge exchange

    if     ( qX ==  2 ) { qM =  1; qB = 1;}
    else if( qX ==  1 ) { qM =  0; qB = 1;}
    else if( qX ==  0 ) { qM =  0; qB = 0;}
    else if( qX == -1 ) { qM = -1; qB = 0;}

    // if     ( qM ==  0 ) pdgM =  pdgM - 100;
    // else if( qM == -1 ) pdgM = -pdgM;

    MesonDecay( lvM, qM); // pdgM ); //

    // else              
    ClusterDecay( lvB, qB ); // continue
  }
}

////////////////////////////////////////////////////////////
//
// Excited barion decay to meson and barion, 
// mass distributions and charge exchange are free parameters 

void G4NeutrinoNucleusModel::MesonDecay( G4LorentzVector & lvX, G4int qX)
{
  G4bool   finB = false;
  G4int    pdgM(0), pdgB(0), i(0), qM(0), qB(0);
  G4double mM(0.), mB(0.), eM(0.), eB(0.), pM(0.), pB(0.);
  G4double mm1(0.), mm22(0.), M1(0.), M2(0.), mX(0.);

  mX = lvX.m();

  G4double mPi = G4ParticleTable::GetParticleTable()->FindParticle(211)->GetPDGMass();

  G4double deltaMr[4] =  { 0.*MeV, 0.*MeV, 100.*MeV, 0.*MeV}; 

  G4ThreeVector   dir(0.,0.,0.);
  G4ThreeVector   bst(0.,0.,0.);
  G4LorentzVector lvM(0.,0.,0.,0.);
  G4LorentzVector lvB(0.,0.,0.,0.);

  for( i = 0; i < fClustNumber; ++i) // check resonance
  {
    if( mX >= fMesMass[i] )
    {
        pdgB = fMesPDG[i];
        // mB = G4ParticleTable::GetParticleTable()->FindParticle(pdgB)->GetPDGMass();
        break;
    }
  }
  if( i == fClustNumber ) // || i == fClustNumber-1 ) // low mass, p || n
  {
    if     ( qX ==  1) { pdgB = 211;  qB =  1;} // pi+   
    else if( qX == 0 ) { pdgB = 111;  qB =  0;} // pi0  
    else if( qX == -1) { pdgB = -211; qB = -1;} // pi-

    return FinalMeson( lvX, qB, pdgB);     
  }
  else if( mX < fMesMass[i] + deltaMr[i] ) // || mX < mPi + mPi ) //
  {
    finB = true; // final barion -> out
    pdgB = fMesPDG[i];
    
    // if     ( qX == 1  && pdgB != 2212)  pdgB = pdgB - 10;
    
    if( qX == 0 )        pdgB = pdgB - 100;
    else if( qX == -1 )  pdgB = -pdgB;

    if( finB ) return FinalMeson( lvX, qX, pdgB ); // out
  }
  // no resonance, try 1->2 decay in COM frame

  // try meson

  mm1 = mPi + 1.*MeV;     // pi+
  mm22 = mX - mPi - 1.*MeV; // mX-n

  if( mm22 <= mm1 ) // out
  {
    if     ( qX ==  1) { pdgB = 211; qB = 1;} // pi+   
    else if( qX == 0 ) { pdgB = 111; qB = 0;} // pi0  
    else if( qX == -1) { pdgB = -211; qB = -1;} // pi-

    return FinalMeson(lvX, qB, pdgB);
  }
  else // try decay -> pion + meson(cluster)
  {
    // G4double sigmaM = 50.*MeV; // 100.*MeV; // 200.*MeV; // 400.*MeV; // 800.*MeV; //
    G4double rand   = G4UniformRand();

    if     ( qX ==  1 ) { qM =  1; qB = 0;}
    else if( qX ==  0 ) { qM =  -1; qB = 1;} // { qM =  0; qB = 0;} //
    else if( qX == -1 ) { qM = -1; qB = 0;}
    /*   
    mM = mPi;
    if(qM == 0) mM = G4ParticleTable::GetParticleTable()->FindParticle(111)->GetPDGMass(); //pi0
    pdgM = fMesPDG[fClustNumber-1];
    */
    // mm1*mm22/( mm1 + rand*(mm22 - mm1) );
    // mM = mm1*mm22/sqrt( mm1*mm1 + rand*(mm22*mm22 - mm1*mm1) );
    // mM = -sigmaM*log( (1.- rand)*exp(-mm22/sigmaM) + rand*exp(-mm1/sigmaM) );
    mM = mm1 + rand*(mm22-mm1);
    // mM = mm1 + 0.9*(mm22-mm1);

    
    for( i = 0; i < fClustNumber; ++i)
    {
      if( mM >= fMesMass[i] )
      {
        pdgM = fMesPDG[i];
        // mM = G4ParticleTable::GetParticleTable()->FindParticle(pdgM)->GetPDGMass();
        break;
      }
    }
    if( i == fClustNumber || i == fClustNumber-1 ) // low mass, p || n
    {
      if     ( qX ==  1) { pdgB = 211; qB = 1;} // pi+   
      else if( qX == 0 ) { pdgB = 111; qB = 0;} // pi0  
      else if( qX == -1) { pdgB = -211; qB = -1;} // pi-

      return FinalMeson( lvX, qB, pdgB);     
    }
    else if( mX < fMesMass[i] + deltaMr[i] ) // || mX < mPi + mPi ) //
    {
      finB = true; // final barion -> out
      pdgB = fMesPDG[i];
    
    // if     ( qX == 1  && pdgB != 2212)  pdgB = pdgB - 10;
    
      if( qX == 0 )        pdgB = pdgB - 100;
      else if( qX == -1 )  pdgB = -pdgB;

      if( finB ) return FinalMeson( lvX, qX, pdgB ); // out
    }
    
    M1 = G4ParticleTable::GetParticleTable()->FindParticle(211)->GetPDGMass()+2.*MeV; // n
    M2 = mX - mM;

    if( M2 <= M1 ) // 
    {
      if     ( qX ==  1) { pdgB = 211; qB = 1;} // pi+   
      else if( qX == 0 ) { pdgB = 111; qB = 0;} // pi0  
      else if( qX == -1) { pdgB = -211; qB = -1;} // pi-

      return FinalMeson(lvX, qB, pdgB);     
    }
    mB = M1 + G4UniformRand()*(M2-M1);
    // mB = -sigmaM*log( (1.- rand)*exp(-M2/sigmaM) + rand*exp(-M1/sigmaM) );
    // mB = M1 + 0.9*(M2-M1);

    dir = G4RandomDirection();
    bst = lvX.boostVector();

    eM  = 0.5*(mX*mX + mM*mM - mB*mB)/mX;
    pM  = sqrt(eM*eM - mM*mM);
    lvM = G4LorentzVector( pM*dir, eM);
    lvM.boost(bst);

    eB  = 0.5*(mX*mX + mB*mB - mM*mM)/mX;
    pB  = sqrt(eB*eB - mB*mB);
    lvB = G4LorentzVector(-pB*dir, eB);
    lvB.boost(bst);

    // G4cout<<mM<<"/"<<mB<<", ";

    // charge exchange

    // if     ( qX ==  2 ) { qM =  1; qB = 1;}

    if     ( qM ==  0 ) pdgM =  pdgM - 100;
    else if( qM == -1 ) pdgM = -pdgM;

    MesonDecay( lvM, qM ); //

    MesonDecay( lvB, qB ); // continue
  }
}

///////////////////////////////////////////////////////////////////////
//
// return final momentum x in the reaction lvX + mI -> mF + mP with momenta p-x, x

G4double G4NeutrinoNucleusModel::FinalMomentum(G4double mI, G4double mF, G4double mP, G4LorentzVector lvX)
{
  G4double result(0.), delta(0.);
  //  G4double mI2 = mI*mI;
  G4double mF2 = mF*mF;
  G4double mP2 = mP*mP;
  G4double eX = lvX.e();
  // G4double mX = lvX.m();
  G4double pX = lvX.vect().mag();
  G4double pX2 = pX*pX;
  G4double sI = eX + mI;
  G4double sI2 = sI*sI;
  G4double B = sI2 - mF2 -pX2 + mP2;
  G4double B2 = B*B;
  G4double a = 4.*(sI2-pX2);
  G4double b = -4.*B*pX;
  G4double c = 4.*sI2*mP2 - B2;
  G4double delta2 = b*b -4.*a*c;

  if( delta2 >= 0. ) delta = sqrt(delta2);

  result = 0.5*(-b-delta)/a;
  // result = 0.5*(-b+delta)/a;

  return result; 
}

/////////////////////////////////////////////////////////////////
//
//

G4double G4NeutrinoNucleusModel::FermiMomentum( G4Nucleus & targetNucleus)
{
  G4int Z  = targetNucleus.GetZ_asInt();
  G4int A  = targetNucleus.GetA_asInt();

  G4double kF(250.*MeV);
  G4double kp = 365.*MeV;
  G4double kn = 231.*MeV;
  G4double t1 = 0.479;
  G4double t2 = 0.526;
  G4double ZpA  = G4double(Z)/G4double(A);
  G4double NpA = 1. - ZpA;

  if      ( Z == 1  && A == 1   ) { kF = 0.;       } // hydrogen ???
  else if ( Z == 1  && A == 2   ) { kF = 87.*MeV;  }
  else if ( Z == 2  && A == 3   ) { kF = 134.*MeV; }
  else if ( Z == 6  && A == 12  ) { kF = 221.*MeV; }
  else if ( Z == 14 && A == 28  ) { kF = 239.*MeV; }
  else if ( Z == 26 && A == 56  ) { kF = 257.*MeV; }
  else if ( Z == 82 && A == 208 ) { kF = 265.*MeV; }
  else 
  {
    kF = kp*ZpA*( 1 - pow( G4double(A), -t1 ) ) + kn*NpA*( 1 - pow( G4double(A), -t2 ) );
  }
  return kF;  
}

/////////////////////////////////////////////////////////////////
//
// sample nucleon momentum of Fermi motion for 1p1h and 2p2h modes

G4double G4NeutrinoNucleusModel::NucleonMomentum( G4Nucleus & targetNucleus)
{
  G4int A     = targetNucleus.GetA_asInt();
  G4double kF = FermiMomentum( targetNucleus);
  G4double mom(0.), kCut = 0.5*GeV; // kCut = 1.*GeV; // kCut = 2.*GeV; // kCut = 4.*GeV; //
  //  G4double cof = 2./GeV;
  // G4double ksi = kF*kF*cof*cof/pi/pi;
  G4double th  = 1.; // 1. - 6.*ksi; //

  if( G4UniformRand() < th || A < 3 )  // 1p1h
  {
    mom = kF*pow( G4UniformRand(), 1./3.); 
  }
  else // 2p2h
  {
    mom  = kF*kCut;
    mom /= kCut - G4UniformRand()*(kCut - kF);
    f2p2h = true;
  }
  return mom;
}

///////////////////////////////////// experimental arrays and get functions ////////////////////////////////////////
//
// Return index of nu/anu energy array corresponding to the neutrino energy

G4int G4NeutrinoNucleusModel::GetEnergyIndex(G4double energy)
{
  G4int i, eIndex = 0;

  for( i = 0; i < fIndex; i++)
  {
    if( energy <= fNuMuEnergy[i]*GeV ) 
    {
      eIndex = i;
      break;
    }
  }
  if( i >= fIndex ) eIndex = fIndex;
  // G4cout<<"eIndex = "<<eIndex<<G4endl;
  return eIndex;
}

/////////////////////////////////////////////////////
//
// nu_mu QE/Tot ratio for index-1, index linear over energy

G4double G4NeutrinoNucleusModel::GetNuMuQeTotRat(G4int index, G4double energy)
{
  G4double ratio(0.);
  // GetMinNuMuEnergy()
  if( index <= 0 || energy < fNuMuEnergy[0] ) ratio = 0.;
  else if (index >= fIndex) ratio = fNuMuQeTotRat[fIndex-1]*fOnePionEnergy[fIndex-1]*GeV/energy;
  else
  {
    G4double x1 = fNuMuEnergy[index-1]*GeV;
    G4double x2 = fNuMuEnergy[index]*GeV;
    G4double y1 = fNuMuQeTotRat[index-1];
    G4double y2 = fNuMuQeTotRat[index];

    if(x1 >= x2) return fNuMuQeTotRat[index];
    else
    {
      G4double angle = (y2-y1)/(x2-x1);
      ratio = y1 + (energy-x1)*angle;
    }
  }
  return ratio;
}

////////////////////////////////////////////////////////

const G4double G4NeutrinoNucleusModel::fNuMuEnergy[50] = 
{
  0.112103, 0.117359, 0.123119, 0.129443, 0.136404, 
  0.144084, 0.152576, 0.161991, 0.172458, 0.184126, 
  0.197171, 0.211801, 0.228261, 0.24684, 0.267887, 
  0.291816, 0.319125, 0.350417, 0.386422, 0.428032, 
  0.47634, 0.532692, 0.598756, 0.676612, 0.768868, 
  0.878812, 1.01062, 1.16963, 1.36271, 1.59876, 
  1.88943, 2.25002, 2.70086, 3.26916, 3.99166, 
  4.91843, 6.11836, 7.6872, 9.75942, 12.5259, 
  16.2605, 21.3615, 28.4141, 38.2903, 52.3062, 
  72.4763, 101.93, 145.6, 211.39, 312.172 
};

////////////////////////////////////////////////////////

const G4double G4NeutrinoNucleusModel::fNuMuQeTotRat[50] = 
{
  // 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 
  // 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 
  // 1., 1., 1., 0.982311, 
  0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98,
  0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98,
  0.97, 0.96, 0.95, 0.93,
  0.917794, 0.850239, 0.780412, 0.709339, 0.638134, 0.568165, 
  0.500236, 0.435528, 0.375015, 0.319157, 0.268463, 0.2232, 0.183284, 
  0.148627, 0.119008, 0.0940699, 0.0733255, 0.0563819, 0.0427312, 0.0319274, 
  0.0235026, 0.0170486, 0.0122149, 0.00857825, 0.00594018, 0.00405037 
};

/////////////////////////////////////////////////////
//
// Return index of one pion array corresponding to the neutrino energy

G4int G4NeutrinoNucleusModel::GetOnePionIndex(G4double energy)
{
  G4int i, eIndex = 0;

  for( i = 0; i < fOnePionIndex; i++)
  {
    if( energy <= fOnePionEnergy[i]*GeV ) 
    {
      eIndex = i;
      break;
    }
  }
  if( i >= fOnePionIndex ) eIndex = fOnePionIndex;
  // G4cout<<"eIndex = "<<eIndex<<G4endl;
  return eIndex;
}

/////////////////////////////////////////////////////
//
// nu_mu 1pi/Tot ratio for index-1, index linear over energy

G4double G4NeutrinoNucleusModel::GetNuMuOnePionProb(G4int index, G4double energy)
{
  G4double ratio(0.);

  if(       index <= 0 || energy < fOnePionEnergy[0] ) ratio = 0.;
  else if ( index >= fOnePionIndex )      ratio = fOnePionProb[fOnePionIndex-1]*fOnePionEnergy[fOnePionIndex-1]*GeV/energy;
  else
  {
    G4double x1 = fOnePionEnergy[index-1]*GeV;
    G4double x2 = fOnePionEnergy[index]*GeV;
    G4double y1 = fOnePionProb[index-1];
    G4double y2 = fOnePionProb[index];

    if( x1 >= x2) return fOnePionProb[index];
    else
    {
      G4double angle = (y2-y1)/(x2-x1);
      ratio = y1 + (energy-x1)*angle;
    }
  }
  return ratio;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

const G4double G4NeutrinoNucleusModel::fOnePionEnergy[58] = 
{

  0.275314, 0.293652, 0.31729, 0.33409, 0.351746, 0.365629, 0.380041, 0.400165, 0.437941, 0.479237, 
  0.504391, 0.537803, 0.588487, 0.627532, 0.686839, 0.791905, 0.878332, 0.987405, 1.08162, 1.16971, 
  1.2982, 1.40393, 1.49854, 1.64168, 1.7524, 1.87058, 2.02273, 2.15894, 2.3654, 2.55792, 2.73017, 
  3.03005, 3.40733, 3.88128, 4.53725, 5.16786, 5.73439, 6.53106, 7.43879, 8.36214, 9.39965, 10.296, 
  11.5735, 13.1801, 15.2052, 17.5414, 19.7178, 22.7462, 25.9026, 29.4955, 33.5867, 39.2516, 46.4716, 
  53.6065, 63.4668, 73.2147, 85.5593, 99.9854 
};


////////////////////////////////////////////////////////////////////////////////////////////////////

const G4double G4NeutrinoNucleusModel::fOnePionProb[58] = 
{
  0.0019357, 0.0189361, 0.0378722, 0.0502758, 0.0662559, 0.0754581, 0.0865008, 0.0987275, 0.124112, 
  0.153787, 0.18308, 0.213996, 0.245358, 0.274425, 0.301536, 0.326612, 0.338208, 0.337806, 0.335948, 
  0.328092, 0.313557, 0.304965, 0.292169, 0.28481, 0.269474, 0.254138, 0.247499, 0.236249, 0.221654, 
  0.205492, 0.198781, 0.182216, 0.162251, 0.142878, 0.128631, 0.116001, 0.108435, 0.0974843, 0.082092, 
  0.0755204, 0.0703121, 0.0607066, 0.0554278, 0.0480401, 0.0427023, 0.0377123, 0.0323248, 0.0298584, 
  0.0244296, 0.0218526, 0.019121, 0.016477, 0.0137309, 0.0137963, 0.0110371, 0.00834028, 0.00686127, 0.00538226
};
 
//
//
///////////////////////////
