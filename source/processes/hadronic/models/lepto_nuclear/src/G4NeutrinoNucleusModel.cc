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
  fEindex = fXindex = fQindex = 0;
  fOnePionIndex = 58;
  fIndex = 50;
  fCascade = fString = fProton = f2p2h = fBreak = false;

  fNuEnergy  = fQ2  = fQtransfer = fXsample = fDp = fTr = 0.;
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

  theMuonMinus = G4MuonMinus::MuonMinus();
  theMuonPlus = G4MuonPlus::MuonPlus();

  // PDG2016: sin^2 theta Weinberg

  fSin2tW = 0.23129; // 0.2312;

  fCutEnergy = 0.; // default value


  /* 
  // G4VPreCompoundModel* ptr;
  // reuse existing pre-compound model as in binary cascade
  
   fPrecoInterface = new  G4GeneratorPrecompoundInterface ;  
 
    if( !fPreCompound )
    {
      G4HadronicInteraction* p =
        G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
      G4VPreCompoundModel* fPreCompound = static_cast<G4VPreCompoundModel*>(p);
      
      if(!fPreCompound)
      {
	fPreCompound = new G4PreCompoundModel();
      }
      fPrecoInterface->SetDeExcitation(fPreCompound);
    }
    fDeExcitation = GetDeExcitation()->GetExcitationHandler();
  */ 
    
  fDeExcitation   = new G4ExcitationHandler();
  fPreCompound    = new G4PreCompoundModel(fDeExcitation);
  fPrecoInterface = new  G4GeneratorPrecompoundInterface ;  
  fPrecoInterface->SetDeExcitation(fPreCompound);
  
  fPDGencoding = 0; // unphysical as default
  fRecoil = nullptr;
     
}


G4NeutrinoNucleusModel::~G4NeutrinoNucleusModel()
{
 if(fPrecoInterface) delete fPrecoInterface;
}


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
  // G4bool FiNucleon(false);
  
  // if      ( qB ==  1 ) pdg = pdgB - 10;
  // else if ( qB ==  0 ) pdg = pdgB - 110;
  // else if ( qB == -1 ) pdg = pdgB - 1110;

  if( pdg == 2212 || pdg == 2112)
  {
    fMr = G4ParticleTable::GetParticleTable()->FindParticle(pdg)->GetPDGMass();
    // FiNucleon = true;
  }
  else fMr = lvB.m();
  
  G4ThreeVector bst = fLVt.boostVector();
  lvB.boost(-bst); // in fLVt rest system
  
  G4double eX = lvB.e();
  G4double det(0.), det2(0.), rM(0.), mX = lvB.m();
  G4ThreeVector dX = (lvB.vect()).unit();
  G4double pX = sqrt(eX*eX-mX*mX);

  if( fRecoil )
  {
    Z  = fRecoil->GetZ_asInt();
    A  = fRecoil->GetA_asInt();
    rM = fRecoil->AtomicMass(A,Z); //->AtomicMass(); // 
    rM = fLVt.m();
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
  det2       = b*b-4.*a*c;
  if( det2 <= 0. ) det = 0.;
  else             det = sqrt(det2);
  G4double dP = 0.5*(-b - det )/a;

  fDp = dP;

  pX -= dP;

  if(pX < 0.) pX = 0.;

  // if( A == 0 ) G4cout<<pX/MeV<<", ";

  eX = sqrt( pX*pX + fMr*fMr );
  G4LorentzVector lvN( pX*dX, eX );
  lvN.boost(bst); // back to lab

  if( pdg == 2212 || pdg == 2112) // nucleons mX >= fMr, dP >= 0
  {
    G4ParticleDefinition* pd2 = G4ParticleTable::GetParticleTable()->FindParticle(pdg); 
    G4DynamicParticle*    dp2 = new G4DynamicParticle( pd2, lvN);
    theParticleChange.AddSecondary( dp2 );

  }
  else // delta resonances
  {
    G4ParticleDefinition* rePart = G4ParticleTable::GetParticleTable()->FindParticle(pdg); 
    G4KineticTrack ddkt( rePart, 0., G4ThreeVector( 0., 0., 0.), lvN);
    G4KineticTrackVector* ddktv = ddkt.Decay();

    G4DecayKineticTracks decay( ddktv );

    for( unsigned int i = 0; i < ddktv->size(); i++ ) // add products to partchange
    {
      G4DynamicParticle * aNew = 
      new G4DynamicParticle( ddktv->operator[](i)->GetDefinition(),
                             ddktv->operator[](i)->Get4Momentum()  );

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
  // dP += G4UniformRand()*10.*MeV;
  G4LorentzVector rec4v(vRecoil, 0.);
  rec4v.boost(bst); // back to lab
  fLVt += rec4v;
  const G4LorentzVector lvTarg = fLVt; // (vRecoil, eRecoil);


  if( fRecoil ) // proton*?
  {
    G4double grM = G4NucleiProperties::GetNuclearMass(A,Z);
    G4double exE = fLVt.m() - grM;
    if( exE < 5.*MeV ) exE = 5.*MeV + G4UniformRand()*10.*MeV;
    
    const G4LorentzVector in4v( G4ThreeVector( 0., 0., 0.), grM );   
    G4Fragment fragment( A, Z, in4v); // lvTarg );
    fragment.SetNumberOfHoles(1);
    fragment.SetExcEnergyAndMomentum( exE, lvTarg );
    
    RecoilDeexcitation(fragment);
  }
  else // momentum?
  {
    theParticleChange.SetLocalEnergyDeposit( fTr ); // eRecoil );
  }
}


///////////////////////////////////////////////////////

void G4NeutrinoNucleusModel::RecoilDeexcitation( G4Fragment& fragment)
{
  G4ReactionProductVector* products = fPreCompound->DeExcite(fragment);

  if( products != NULL )
  {
    G4ReactionProductVector::iterator iter;

    for( iter = products->begin(); iter != products->end(); ++iter )
    {
      G4DynamicParticle * aNewDP =
                    new G4DynamicParticle((*iter)->GetDefinition(),
                            (*iter)->GetTotalEnergy(),
                            (*iter)->GetMomentum());
      /*      
      // G4HadSecondary aNew = G4HadSecondary(aNewDP);

      G4double time=(*iter)->GetFormationTime();

      if( time < 0.0) { time = 0.0; }

      aNew.SetTime(time);// (timePrimary + time);
      aNew.SetCreatorModelType((*iter)->GetCreatorModel());
*/
      // G4cout<<aNewDP->GetDefinition()->GetParticleName()<<", "<<aNewDP->Get4Momentum()<<G4endl;

      theParticleChange.AddSecondary(aNewDP);
    }
  }
  return;
}


///////////////////////////////////////////
//
// Fragmentation of lvX directly to pion and recoil nucleus (A,Z): fLVh + fLVt -> pi + A*

void G4NeutrinoNucleusModel::CoherentPion( G4LorentzVector & lvP, G4int pdgP, G4Nucleus & targetNucleus)
{
  G4int A(0), Z(0), pdg = pdgP;
  fLVcpi = G4LorentzVector(0.,0.,0.,0.);

  G4double  rM(0.), mN(938.), mI(0.), det(0.), det2(0.); 

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
    bst = fLVt.boostVector();  // to recoil rest frame
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
  det2 = b*b-4.*a*c;
  if(det2 > 0.) det = sqrt(det2);
  G4double dP = 0.5*(-b - det )/a;

  dP = FinalMomentum( mI, rM, fMr, lvP);

  // G4cout<<dP<<", ";
  pX -= dP;
  if( pX < 0. ) pX = 0.;
  
  eX = sqrt( dP*dP + fMr*fMr );
  G4LorentzVector lvN( dP*dX, eX );

  if( A > 1 ) lvN.boost(bst); // back to lab

  fLVcpi = lvN;

  G4ParticleDefinition* pd2 = G4ParticleTable::GetParticleTable()->FindParticle(pdg); 
  G4DynamicParticle*    dp2 = new G4DynamicParticle( pd2, lvN);
  theParticleChange.AddSecondary( dp2 ); // coherent pion

  // recoil nucleus

  G4double eRecoil = sqrt( rM*rM + pX*pX );
  G4ThreeVector vRecoil(pX*dX);
  G4LorentzVector lvTarg1( vRecoil, eRecoil);
  lvTarg1.boost(bst);
  
  const G4LorentzVector lvTarg = lvTarg1;

  if( A > 1 ) // recoil target nucleus*
  {
    G4double grM = G4NucleiProperties::GetNuclearMass(A,Z);
    G4double exE = fLVt.m() - grM;
    
    if( exE < 5.*MeV ) exE = 5.*MeV + G4UniformRand()*10.*MeV; // vmg???
    
    const G4LorentzVector in4v( G4ThreeVector( 0., 0., 0.), grM );   
    G4Fragment fragment( A, Z, in4v); // lvTarg );
    fragment.SetNumberOfHoles(1);
    fragment.SetExcEnergyAndMomentum( exE, lvTarg );
    
    RecoilDeexcitation(fragment);
  }
  else // recoil target proton 
  {
    G4double eTkin = eRecoil - rM;
    G4double eTh   = 0.01*MeV; // 10.*MeV;
    
    if( eTkin > eTh )
    {
      G4DynamicParticle * aSec = new G4DynamicParticle( G4Proton::Proton(), lvTarg);
      theParticleChange.AddSecondary(aSec);
    }
    else theParticleChange.SetLocalEnergyDeposit( eTkin );
  }
  return;
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


//////////////////////////////////////////////////////
//
// Excitation energy according Bodek

G4double G4NeutrinoNucleusModel::GetEx( G4int A, G4bool fP )
{
  G4double eX(10.*MeV), a1(0.), a2(0.), e1(0.), e2(0.), aa = G4double(A);
  G4int i(0);
  const G4int maxBin = 12;
  
  G4double refA[maxBin] = { 2., 6., 12., 16., 27., 28., 40., 50., 56., 58., 197., 208. };

  G4double pEx[maxBin] = { 0., 12.2, 10.1, 10.9, 21.6, 12.4, 17.8, 17., 19., 16.8, 19.5, 14.7 };  

  G4double nEx[maxBin] = { 0., 12.2, 10., 10.2, 21.6, 12.4, 21.8, 17., 19., 16.8, 19.5, 16.9 };

  G4DataVector dE(12,0.);

  if(fP) for( i = 0; i < maxBin; ++i ) dE[i] = pEx[i];
  else                                 dE[i] = nEx[i];

  for( i = 0; i < maxBin; ++i )
  {
    if( aa <= refA[i] ) break;
  }
  if( i >= maxBin ) eX = dE[maxBin-1];
  else if( i <= 0 ) eX = dE[0];
  else
  {
    a1 = refA[i-1];
    a2 = refA[i];
    e1 = dE[i-1];
    e2 = dE[i];
    if (a1 == a2 || e1 == e2 ) eX = dE[i];
    else eX = e1 + (e2-e1)*(aa-a1)/(a2-a1);
  }
  return eX;
}



///////////////////////////////////////////////////////
//
// Two G-function sampling for the nucleon momentum  

G4double G4NeutrinoNucleusModel::GgSampleNM(G4Nucleus & nucl)
{
  f2p2h = false;
  G4double /* distr(0.), tail(0.), */ shift(1.), xx(1.), mom(0.), th(0.1);
  G4double kF = FermiMomentum( nucl);
  G4double momMax = 2.*kF; //  1.*GeV; //  1.*GeV; // 
  G4double aa = 5.5;
  G4double ll = 6.0; //  6.5; //

  G4int A = nucl.GetA_asInt();

  if( A <= 12) th = 0.1;
  else
  {
    // th = 0.1/(1.+log(G4double(A)/12.));
    th = 1.2/( G4double(A) + 1.35*log(G4double(A)/12.) );
  }
  shift = 0.99; // 0.95; //
  xx = mom/shift/kF;

  G4double rr = G4UniformRand();

  if( rr > th )
  {
    aa = 5.5;
    
    if( A <= 12 ) ll = 6.0;
    else
    {
      ll = 6.0 + 1.35*log(G4double(A)/12.);
    }
    xx = RandGamma::shoot(aa,ll);
    shift = 0.99;
    mom = xx*shift*kF;
  }
  else
  {
    f2p2h = true;
    aa = 6.5;
    ll = 6.5;
    xx = RandGamma::shoot(aa,ll);
    shift = 2.5;
    mom = xx*shift*kF;
  }
  if( mom > momMax ) mom = G4UniformRand()*momMax;
  if( mom > 2.*kF  ) f2p2h = true;

  // mom = 0.;

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
