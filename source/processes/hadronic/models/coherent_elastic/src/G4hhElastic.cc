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
//
// Physics model class G4hhElastic 
//
//
// G4 Model: qQ hadron hadron elastic scattering with 4-momentum balance
//                         
// 02.05.2014 V. Grichine 1-st version
//

#include "G4hhElastic.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4NucleiProperties.hh"

#include "Randomize.hh"
#include "G4Integrator.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"

#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsFreeVector.hh"

#include "G4HadronNucleonXsc.hh"

#include "G4Pow.hh"

#include "G4HadronicParameters.hh"

using namespace std;


/////////////////////////////////////////////////////////////////////////
//
// Tracking constructor. Target is proton


G4hhElastic::G4hhElastic() 
  : G4HadronElastic("HadrHadrElastic")
{
  SetMinEnergy( 1.*GeV );
  SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  verboseLevel = 0;
  lowEnergyRecoilLimit = 100.*keV;  
  lowEnergyLimitQ  = 0.0*GeV;  
  lowEnergyLimitHE = 0.0*GeV;  
  lowestEnergyLimit= 0.0*keV;  
  plabLowLimit     = 20.0*MeV;

  fRhoReIm=fSigmaTot=fOptRatio=fSpp=fPcms=0.0;
  fInTkin=0;
  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  thePionPlus = G4PionPlus::PionPlus();
  thePionMinus= G4PionMinus::PionMinus();

  fTarget  = G4Proton::Proton();
  fProjectile  = 0;
  fHadrNuclXsc = new G4HadronNucleonXsc();

  fEnergyBin = 200;
  fBinT  = 514; // 514; // 500; // 200;

  fEnergyVector =  new G4PhysicsLogVector( theMinEnergy, theMaxEnergy, fEnergyBin );

  fTableT = 0;
  fOldTkin = 0.;
  SetParameters();

  Initialise();
}


/////////////////////////////////////////////////////////////////////////
//
// test constructor


G4hhElastic::G4hhElastic( G4ParticleDefinition* target, G4ParticleDefinition* projectile, G4double plab) 
  : G4HadronElastic("HadrHadrElastic")
{
  SetMinEnergy( 1.*GeV );
  SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  verboseLevel         = 0;
  lowEnergyRecoilLimit = 100.*keV;  
  lowEnergyLimitQ      = 0.0*GeV;  
  lowEnergyLimitHE     = 0.0*GeV;  
  lowestEnergyLimit    = 0.0*keV;  
  plabLowLimit         = 20.0*MeV;

  fRhoReIm=fSigmaTot=fOptRatio=fSpp=fPcms=0.0;
  fInTkin=0;
  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  thePionPlus = G4PionPlus::PionPlus();
  thePionMinus= G4PionMinus::PionMinus();

  fTarget      = target;
  fProjectile  = projectile;
  fMassTarg   = fTarget->GetPDGMass();
  fMassProj   = fProjectile->GetPDGMass();
  fMassSum2   = (fMassTarg+fMassProj)*(fMassTarg+fMassProj);
  fMassDif2   = (fMassTarg-fMassProj)*(fMassTarg-fMassProj);
  fHadrNuclXsc = new G4HadronNucleonXsc();

  fEnergyBin = 200;
  fBinT      = 514; // 200;

  fEnergyVector =  new G4PhysicsLogVector( theMinEnergy, theMaxEnergy, fEnergyBin );
  fTableT       = 0;
  fOldTkin = 0.;


  SetParameters();
  SetParametersCMS( plab);
}


/////////////////////////////////////////////////////////////////////////
//
// constructor used for low mass diffraction


G4hhElastic::G4hhElastic( G4ParticleDefinition* target, G4ParticleDefinition* projectile) 
  : G4HadronElastic("HadrHadrElastic")
{
  SetMinEnergy( 1.*GeV );
  SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  verboseLevel = 0;
  lowEnergyRecoilLimit = 100.*keV;  
  lowEnergyLimitQ  = 0.0*GeV;  
  lowEnergyLimitHE = 0.0*GeV;  
  lowestEnergyLimit= 0.0*keV;  
  plabLowLimit     = 20.0*MeV;

  fRhoReIm=fSigmaTot=fOptRatio=fSpp=fPcms=0.0;
  fInTkin=0;

  fTarget = target; // later vmg
  fProjectile = projectile;
  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  thePionPlus = G4PionPlus::PionPlus();
  thePionMinus= G4PionMinus::PionMinus();

  fTarget  = G4Proton::Proton(); // later vmg
  fMassTarg   = fTarget->GetPDGMass();
  fMassProj   = fProjectile->GetPDGMass();
  fMassSum2   = (fMassTarg+fMassProj)*(fMassTarg+fMassProj);
  fMassDif2   = (fMassTarg-fMassProj)*(fMassTarg-fMassProj);
  fHadrNuclXsc = new G4HadronNucleonXsc();

  fEnergyBin = 200;
  fBinT  = 514; // 514; // 500; // 200;

  fEnergyVector =  new G4PhysicsLogVector( theMinEnergy, theMaxEnergy, fEnergyBin );

  fTableT = 0;
  fOldTkin = 0.;

  SetParameters();
}



//////////////////////////////////////////////////////////////////////////////
//
// Destructor

G4hhElastic::~G4hhElastic()
{
  if ( fEnergyVector ) {
    delete fEnergyVector;
    fEnergyVector = 0;
  }

  for ( std::vector<G4PhysicsTable*>::iterator it = fBankT.begin();
        it != fBankT.end(); ++it ) {
    if ( (*it) ) (*it)->clearAndDestroy();
    delete *it;
    *it = 0;
  }
  fTableT = 0;
  if(fHadrNuclXsc) delete fHadrNuclXsc;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////  Table preparation and reading ////////////////////////


//////////////////////////////////////////////////////////////////////////////
//
// Initialisation for given particle on the proton target

void G4hhElastic::Initialise() 
{
  // pp,pn

  fProjectile = G4Proton::Proton();
  BuildTableT(fTarget, fProjectile);
  fBankT.push_back(fTableT); // 0

  // pi+-p
  
  fProjectile = G4PionPlus::PionPlus();
  BuildTableT(fTarget, fProjectile);
  fBankT.push_back(fTableT);  // 1
  //K+-p
  fProjectile = G4KaonPlus::KaonPlus();
  BuildTableT(fTarget, fProjectile);
  fBankT.push_back(fTableT);  // 2
  
}

///////////////////////////////////////////////////////////////////////////////
//
// Build for given particle and proton table of momentum transfers.

void G4hhElastic::BuildTableT( G4ParticleDefinition* target, G4ParticleDefinition* projectile) // , G4double plab) 
{
  G4int iTkin, jTransfer;
  G4double plab, Tkin, tMax;
  G4double t1, t2, dt, delta = 0., sum = 0.;

  fTarget     = target;
  fProjectile = projectile;
  fMassTarg   = fTarget->GetPDGMass();
  fMassProj   = fProjectile->GetPDGMass();
  fMassSum2   = (fMassTarg+fMassProj)*(fMassTarg+fMassProj);
  fMassDif2   = (fMassTarg-fMassProj)*(fMassTarg-fMassProj);

  G4Integrator<G4hhElastic,G4double(G4hhElastic::*)(G4double)> integral;
  // G4HadronNucleonXsc*          hnXsc = new G4HadronNucleonXsc();
  fTableT = new G4PhysicsTable(fEnergyBin);

  for( iTkin = 0; iTkin < fEnergyBin; iTkin++)
  {
    Tkin  = fEnergyVector->GetLowEdgeEnergy(iTkin);
    plab  = std::sqrt( Tkin*( Tkin + 2*fMassProj ) );
    // G4DynamicParticle*  theDynamicParticle = new G4DynamicParticle(projectile,
    //                                           G4ParticleMomentum(0.,0.,1.),
    //                                           Tkin);
    // fSigmaTot = fHadrNuclXsc->GetHadronNucleonXscNS( theDynamicParticle, target );

    SetParametersCMS( plab );

    tMax  = 4.*fPcms*fPcms;
    if( tMax > 15.*GeV*GeV ) tMax = 15.*GeV*GeV; // Check vs. energy ???

    G4PhysicsFreeVector* vectorT = new G4PhysicsFreeVector(fBinT-1);
    sum = 0.;
    dt  = tMax/fBinT;

    // for(j = 1; j < fBinT; j++)

    for( jTransfer = fBinT-1; jTransfer >= 1; jTransfer--)
    {
      t1 = dt*(jTransfer-1);
      t2 = t1 + dt;

      if( fMassProj > 900.*MeV ) // pp, pn
      {
        delta = integral.Legendre10(this, &G4hhElastic::GetdsdtF123, t1, t2);
        // delta = integral.Legendre96(this, &G4hhElastic::GetdsdtF123, t1, t2);
      }
      else // pi+-p, K+-p
      {
        delta = integral.Legendre10(this, &G4hhElastic::GetdsdtF123qQgG, t1, t2);
        // delta = integral.Legendre96(this, &G4hhElastic::GetdsdtF123qQgG, t1, t2);
      }
      sum += delta;
      vectorT->PutValue( jTransfer-1, t1, sum ); // t2
    }
    // vectorT->PutValue( fBinT-1, dt*(fBinT-1), 0. ); // t2
    fTableT->insertAt( iTkin, vectorT );
    // delete theDynamicParticle;
  }
  // delete hnXsc;

  return;
}

////////////////////////////////////////////////////////////////////////////
//
// Return inv momentum transfer -t > 0 from initialisation table

G4double G4hhElastic::SampleInvariantT( const G4ParticleDefinition* aParticle, G4double p, 
                                               G4int, G4int )
{
  G4int iTkin, iTransfer;
  G4double t, t2, position, m1 = aParticle->GetPDGMass();
  G4double Tkin = std::sqrt(m1*m1+p*p) - m1;

  if( aParticle == G4Proton::Proton() || aParticle == G4Neutron::Neutron() )
  {
    fTableT = fBankT[0];
  }
  if( aParticle == G4PionPlus::PionPlus() || aParticle == G4PionMinus::PionMinus() )
  {
    fTableT = fBankT[1];
  }
  if( aParticle == G4KaonPlus::KaonPlus() || aParticle == G4KaonMinus::KaonMinus() )
  {
    fTableT = fBankT[2];
  }

  G4double delta    = std::abs(Tkin - fOldTkin)/(Tkin + fOldTkin);
  G4double deltaMax = 1.e-2;

  if ( delta < deltaMax ) iTkin = fInTkin; 
  else
  {  
    for( iTkin = 0; iTkin < fEnergyBin; iTkin++)
    {
      if( Tkin < fEnergyVector->GetLowEdgeEnergy(iTkin) ) break;
    }
  }
  if ( iTkin >= fEnergyBin ) iTkin = fEnergyBin-1;   // Tkin is more then theMaxEnergy
  if ( iTkin < 0 )           iTkin = 0; // against negative index, Tkin < theMinEnergy

  fOldTkin = Tkin;
  fInTkin = iTkin;

  if (iTkin == fEnergyBin -1 || iTkin == 0 )   // the table edges
  {
    position = (*(*fTableT)(iTkin))(0)*G4UniformRand();

    // G4cout<<"position = "<<position<<G4endl;

    for(iTransfer = 0; iTransfer < fBinT-1; iTransfer++)
    {
      if( position >= (*(*fTableT)(iTkin))(iTransfer) ) break;
    }
    if (iTransfer >= fBinT-1) iTransfer = fBinT-2;

    // G4cout<<"iTransfer = "<<iTransfer<<G4endl;

    t = GetTransfer(iTkin, iTransfer, position);

    // G4cout<<"t = "<<t<<G4endl;
  }
  else  // Tkin inside between energy table edges
  {
    // position = (*(*fTableT)(iTkin))(fBinT-2)*G4UniformRand();
    position = (*(*fTableT)(iTkin))(0)*G4UniformRand();

    // G4cout<<"position = "<<position<<G4endl;

    for(iTransfer = 0; iTransfer < fBinT-1; iTransfer++)
    {
      // if( position < (*(*fTableT)(iTkin))(iTransfer) ) break;
      if( position >= (*(*fTableT)(iTkin))(iTransfer) ) break;
    }
    if (iTransfer >= fBinT-1) iTransfer = fBinT-2;

    // G4cout<<"iTransfer = "<<iTransfer<<G4endl;

    t2  = GetTransfer(iTkin, iTransfer, position);
    return t2;
    /*
    G4double t1, E1, E2, W, W1, W2;
    // G4cout<<"t2 = "<<t2<<G4endl;

    E2 = fEnergyVector->GetLowEdgeEnergy(iTkin);

    // G4cout<<"E2 = "<<E2<<G4endl;
    
    iTkin--;
    
    // position = (*(*fTableT)(iTkin))(fBinT-2)*G4UniformRand();

    // G4cout<<"position = "<<position<<G4endl;

    for(iTransfer = 0; iTransfer < fBinT-1; iTransfer++)
    {
      // if( position < (*(*fTableT)(iTkin))(iTransfer) ) break;
      if( position >= (*(*fTableT)(iTkin))(iTransfer) ) break;
    }
    if (iTransfer >= fBinT-1) iTransfer = fBinT-2;
    
    t1  = GetTransfer(iTkin, iTransfer, position);

    // G4cout<<"t1 = "<<t1<<G4endl;

    E1 = fEnergyVector->GetLowEdgeEnergy(iTkin);

    // G4cout<<"E1 = "<<E1<<G4endl;

    W  = 1.0/(E2 - E1);
    W1 = (E2 - Tkin)*W;
    W2 = (Tkin - E1)*W;

    t = W1*t1 + W2*t2;
    */
  }
  return t;
}


////////////////////////////////////////////////////////////////////////////
//
// Return inv momentum transfer -t > 0 from initialisation table 

G4double G4hhElastic::SampleBisectionalT( const G4ParticleDefinition* aParticle, G4double p)
{
  G4int iTkin, iTransfer;
  G4double t,  position, m1 = aParticle->GetPDGMass();
  G4double Tkin = std::sqrt(m1*m1+p*p) - m1;

  if( aParticle == G4Proton::Proton() || aParticle == G4Neutron::Neutron() )
  {
    fTableT = fBankT[0];
  }
  if( aParticle == G4PionPlus::PionPlus() || aParticle == G4PionMinus::PionMinus() )
  {
    fTableT = fBankT[1];
  }
  if( aParticle == G4KaonPlus::KaonPlus() || aParticle == G4KaonMinus::KaonMinus() )
  {
    fTableT = fBankT[2];
  }
  G4double delta    = std::abs(Tkin - fOldTkin)/(Tkin + fOldTkin);
  G4double deltaMax = 1.e-2;

  if ( delta < deltaMax ) iTkin = fInTkin; 
  else
  {  
    for( iTkin = 0; iTkin < fEnergyBin; iTkin++ )
    {
      if( Tkin < fEnergyVector->GetLowEdgeEnergy(iTkin) ) break;
    }
  }
  if ( iTkin >= fEnergyBin ) iTkin = fEnergyBin-1;   // Tkin is more then theMaxEnergy
  if ( iTkin < 0 )           iTkin = 0; // against negative index, Tkin < theMinEnergy

  fOldTkin = Tkin;
  fInTkin = iTkin;

  if (iTkin == fEnergyBin -1 || iTkin == 0 )   // the table edges
  {
    position = (*(*fTableT)(iTkin))(0)*G4UniformRand();

    for(iTransfer = 0; iTransfer < fBinT-1; iTransfer++)
    {
      if( position >= (*(*fTableT)(iTkin))(iTransfer) ) break;
    }
    if (iTransfer >= fBinT-1) iTransfer = fBinT-2;

    t = GetTransfer(iTkin, iTransfer, position);


  }
  else  // Tkin inside between energy table edges
  {
    G4double rand = G4UniformRand();
    position = (*(*fTableT)(iTkin))(0)*rand;

    //
    // (*fTableT)(iTkin)->GetLowEdgeEnergy(fBinT-2);
    G4int sTransfer = 0, fTransfer =  fBinT - 2, dTransfer = fTransfer - sTransfer;
    G4double y2;

    for( iTransfer = 0; iTransfer < fBinT - 1; iTransfer++ )
    {
      // dTransfer %= 2;
      dTransfer /= 2;
      // dTransfer *= 0.5;
      y2 = (*(*fTableT)(iTkin))( sTransfer + dTransfer );

      if( y2 > position ) sTransfer += dTransfer;

      // if( dTransfer <= 1 ) break;
      if( dTransfer < 1 ) break;
    }
    t = (*fTableT)(iTkin)->GetLowEdgeEnergy(sTransfer); // +(-0.5+rand)*(*fTableT)(iTkin)->GetLowEdgeEnergy(3);
  }
  return t;
}


///////////////////////////////////////////////////////////////////////////////
//
// Build for given particle and proton table of momentum transfers.

void G4hhElastic::BuildTableTest( G4ParticleDefinition* target, G4ParticleDefinition* projectile, G4double plab) 
{
  G4int jTransfer;
  G4double tMax; // , sQq, sQG;
  G4double t1, t2, dt, delta = 0., sum = 0. ; // , threshold;

  fTarget     = target;
  fProjectile = projectile;
  fMassTarg =  fTarget->GetPDGMass();
  fMassProj =  fProjectile->GetPDGMass();
  fMassSum2 = (fMassTarg+fMassProj)*(fMassTarg+fMassProj);
  fMassDif2 = (fMassTarg-fMassProj)*(fMassTarg-fMassProj);
  fSpp  = fMassProj*fMassProj + fMassTarg*fMassTarg + 2.*fMassTarg*std::sqrt(plab*plab + fMassProj*fMassProj);
  fPcms = std::sqrt( (fSpp - fMassSum2)*(fSpp - fMassDif2)/4./fSpp);

  G4cout<<"fMassTarg = "<<fMassTarg<<" MeV; fMassProj = "<<fMassProj<<" MeV"<<G4endl;
  tMax = 4.*fPcms*fPcms;
  if( tMax > 15.*GeV*GeV ) tMax = 15.*GeV*GeV; // Check vs. energy ???

 
  G4Integrator<G4hhElastic,G4double(G4hhElastic::*)(G4double)> integral;
  fTableT = new G4PhysicsTable(1);
  G4PhysicsFreeVector* vectorT = new G4PhysicsFreeVector(fBinT-1);

  sum = 0.;
  dt = tMax/G4double(fBinT);
  G4cout<<"s = "<<std::sqrt(fSpp)/GeV<<" GeV; fPcms = "<<fPcms/GeV
	<<" GeV; qMax = "<<tMax/GeV/GeV<<" GeV2; dt = "<<dt/GeV/GeV<<" GeV2"<<G4endl;

  // G4cout<<"fRA = "<<fRA*GeV<<"; fRB = "<<fRB*GeV<<G4endl;

    // for(jTransfer = 1; jTransfer < fBinT; jTransfer++)
  for( jTransfer = fBinT-1; jTransfer >= 1; jTransfer-- )
  {
    t1 = dt*(jTransfer-1);
    t2 = t1 + dt;

    if( fMassProj > 900.*MeV ) // pp, pn
    {
      delta = integral.Legendre10(this, &G4hhElastic::GetdsdtF123, t1, t2);
      // threshold = integral.Legendre96(this, &G4hhElastic::GetdsdtF123, t1, tMax);
    }
    else // pi+-p, K+-p
    {
      delta = integral.Legendre10(this, &G4hhElastic::GetdsdtF123qQgG, t1, t2);
      // threshold = integral.Legendre96(this, &G4hhElastic::GetdsdtF123qQgG, t1, tMax);
      // delta = integral.Legendre96(this, &G4hhElastic::GetdsdtF123, t1, t2);
    }
    sum += delta;
    // G4cout<<delta<<"\t"<<sum<<"\t"<<threshold<<G4endl;   
  
    // sQq = GetdsdtF123(q1);
    // sQG = GetdsdtF123qQgG(q1);
    // G4cout<<q1/GeV<<"\t"<<sQG*GeV*GeV/millibarn<<"\t"<<sQq*GeV*GeV/millibarn<<G4endl;
    // G4cout<<"sum = "<<sum<<", "; 
    
    vectorT->PutValue( jTransfer-1, t1, sum );  // t2
  }
  // vectorT->PutValue( fBinT-1, dt*(fBinT-1), 0. ); // t2
  fTableT->insertAt( 0, vectorT );
  fBankT.push_back( fTableT );  // 0

  // for(jTransfer = 0; jTransfer < fBinT-1; jTransfer++) 
  //   G4cout<<(*(*fTableT)(0))(jTransfer)/sum<<"\t\t"<<G4Pow::GetInstance()->powN(2.,-jTransfer)<<G4endl;
  
  return;
}


////////////////////////////////////////////////////////////////////////////
//
// Return inv momentum transfer -t > 0 from initialisation table

G4double G4hhElastic::SampleTest(G4double tMin ) // const G4ParticleDefinition* aParticle,  )
{
  G4int iTkin, iTransfer, iTmin;
  G4double t, position;
  // G4double qMin = std::sqrt(tMin);

  fTableT = fBankT[0];
  iTkin = 0;

  for(iTransfer = 0; iTransfer < fBinT-1; iTransfer++)
  {
    // if( qMin <= (*fTableT)(iTkin)->GetLowEdgeEnergy(iTransfer) ) break;
    if( tMin <= (*fTableT)(iTkin)->GetLowEdgeEnergy(iTransfer) ) break;
  }
  iTmin = iTransfer-1;
  if(iTmin < 0 ) iTmin = 0;

  position = (*(*fTableT)(iTkin))(iTmin)*G4UniformRand();

  for( iTmin = 0; iTransfer < fBinT-1; iTransfer++)
  {
    if( position > (*(*fTableT)(iTkin))(iTransfer) ) break;
  }
  if (iTransfer >= fBinT-1) iTransfer = fBinT-2;

  t = GetTransfer(iTkin, iTransfer, position);

  return t;
}


/////////////////////////////////////////////////////////////////////////////////
//
// Check with PAI sampling

G4double 
G4hhElastic:: GetTransfer( G4int iTkin, G4int iTransfer, G4double position )
{
  G4double x1, x2, y1, y2, randTransfer, delta, mean, epsilon = 1.e-6;

  if( iTransfer == 0 )
  {
    randTransfer = (*fTableT)(iTkin)->GetLowEdgeEnergy(iTransfer);
    // iTransfer++;
  }
  else
  {
    if ( iTransfer >= G4int((*fTableT)(iTkin)->GetVectorLength()) )
    {
      iTransfer = G4int((*fTableT)(iTkin)->GetVectorLength() - 1);
    }
    y1 = (*(*fTableT)(iTkin))(iTransfer-1);
    y2 = (*(*fTableT)(iTkin))(iTransfer);

    x1 = (*fTableT)(iTkin)->GetLowEdgeEnergy(iTransfer-1);
    x2 = (*fTableT)(iTkin)->GetLowEdgeEnergy(iTransfer);

    delta = y2 - y1;
    mean  = y2 + y1;

    if ( x1 == x2 ) randTransfer = x2;
    else
    {
      // if ( y1 == y2 ) 
      if ( delta < epsilon*mean ) 
           randTransfer = x1 + ( x2 - x1 )*G4UniformRand();
      else randTransfer = x1 + ( position - y1 )*( x2 - x1 )/delta; // ( y2 - y1 );
    }
  }
  return randTransfer;
}

const G4double G4hhElastic::theNuclNuclData[19][6] = 
{
  // sqrt(fSpp) in GeV, fRA in 1/GeV, fRB in 1/GeV, fBq, fBQ, fImCof 

  { 2.76754,	4.8,	4.8,	0.05,	0.742441,	10.5 }, // pp 3GeV/c
  { 3.07744,	5.4,	5.4,	0.02,	0.83818,	6.5  }, // pp 4GeV/c
  { 3.36305,	5.2,	5.2,	0.02,	0.838893,	7.5  }, // np 5GeV/c
  { 4.32941,	6,	6,	0.03,	0.769389,	7.5  }, // np 9 GeV/c
  { 4.62126,	6,	6,	0.03,	0.770111,	6.5  }, // pp 10.4 GeV/c

  { 5.47416,	4.5,	4.5,	0.03,	0.813185,	7.5  }, // np 15 GeV/c
  { 6.15088,	6.5,	6.5,	0.02,	0.799539,	6.5  }, // pp 19.2 GeV/c
  { 6.77474,	5.2,	5.2,	0.03,	0.784901,	7.5  }, // np 23.5 GeV/c
  { 9.77775,	7,	7,	0.03,	0.742531,	6.5  }, // pp 50 GeV/c
                // {9.77775,	7,	7,	0.011,	0.84419,	4.5 }, // pp 50 GeV/c
  { 10.4728,	5.2,	5.2,	0.03,	0.780439,	7.5  }, // np 57.5 GeV/c

  { 13.7631,	7,	7,	0.008,	0.8664,	        5.0  }, // pp 100 GeV/c
  { 19.4184,	6.8,	6.8,	0.009,	0.861337,	2.5  }, // pp 200 GeV/c
  { 23.5,	6.8,	6.8,	0.007,	0.878112,	1.5  }, // pp 23.5 GeV
                // {24.1362,	6.4,	6.4,	0.09,	0.576215,	7.5 }, // np 309.5 GeV/c
  { 24.1362,	7.2,	7.2,	0.008,	0.864745,	5.5  },
  { 52.8,	6.8,	6.8,	0.008,	0.871929,	1.5  }, // pp 58.2 GeV

  { 546,        7.4,	7.4,	0.013,	0.845877,	5.5  }, // pb-p 546 GeV
  { 1960,	7.8,	7.8,	0.022,	0.809062,	7.5  }, // pb-p 1960 GeV
  { 7000,	8,	8,	0.024,	0.820441,	5.5  },  // pp TOTEM 7 TeV
  { 13000,	8.5,	8.5,	0.03,	0.796721,	10.5  }  // pp TOTEM 13 TeV
 
};

//////////////////////////////////////////////////////////////////////////////////

const G4double G4hhElastic::thePiKaNuclData[8][6] = 
{
  // sqrt(fSpp) in GeV, fRA in 1/GeV, fRB in 1/GeV, fBq, fBQ, fImCof 

  { 2.5627,     3.8,    3.3,    0.22,   0.222,          1.5 }, // pipp 3.017 GeV/c
  { 2.93928,    4.3,    3.8,    0.2,    0.250601,       1.3 }, // pipp 4.122 GeV/c
  { 3.22326,	4.8,	4.3,	0.13,	0.32751,	2.5 }, // pipp 5.055 GeV/c
  { 7.80704,	5.5,	5,	0.13,	0.340631,	2.5 }, // pipp 32 GeV/c
  { 9.7328,	5,	4.5,	0.05,	0.416319,	5.5  }, // pipp 50 GeV/c

  { 13.7315,	5.3,	4.8,	0.05,	0.418426,	5.5  }, // pipp 100 GeV/c
  { 16.6359,	6.3,	5.8,	0.05,	0.423817,	5.5  }, // pipp 147 GeV/c
  { 19.3961,	5,	4.5,	0.05,	0.413477,	3.5  } // pimp 200 GeV/c
 
};

//
//
/////////////////////////////////////////////////////////////////////////////////
