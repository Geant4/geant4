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
// $Id$
//
// G4LMsdGenerator
//
//

#include "G4DynamicParticle.hh"
#include "G4LMsdGenerator.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4IonTable.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadFinalState.hh"
#include "G4KineticTrack.hh"
#include "G4DecayKineticTracks.hh"
#include "G4KineticTrackVector.hh"
#include "G4Log.hh"


G4LMsdGenerator::G4LMsdGenerator(const G4String& name)
    : G4HadronicInteraction(name)

{
  fPDGencoding = 0;

  // theParticleChange = new G4HadFinalState;
}

G4LMsdGenerator::~G4LMsdGenerator()
{
  // delete theParticleChange;
}

void G4LMsdGenerator::ModelDescription(std::ostream& outFile) const
{
  outFile << GetModelName() <<" consists of a " 
	  << " string model and a stage to de-excite the excited nuclear fragment."
	  << "\n<p>"
	  << "The string model simulates the interaction of\n"
          << "an incident hadron with a nucleus, forming \n"
          << "excited strings, decays these strings into hadrons,\n"
          << "and leaves an excited nucleus. \n"
          << "<p>The string model:\n";
}

/////////////////////////////////////////////////////////////////
//
// Particle and kinematical limitation od diffraction dissociation

G4bool   
G4LMsdGenerator::IsApplicable( const G4HadProjectile& aTrack, 
                                      G4Nucleus& targetNucleus )
{
  G4bool applied = false;

  if( ( aTrack.GetDefinition() == G4Proton::Proton() || 
	aTrack.GetDefinition() == G4Neutron::Neutron()  ) && 
        targetNucleus.GetA_asInt() >= 1 &&
      // aTrack.GetKineticEnergy() > 1800*CLHEP::MeV 
        aTrack.GetKineticEnergy() > 300*CLHEP::MeV 
    ) //  750*CLHEP::MeV )   
  {
    applied = true;
  }
  else if( ( aTrack.GetDefinition() == G4PionPlus::PionPlus() || 
	     aTrack.GetDefinition() == G4PionMinus::PionMinus()  ) && 
             targetNucleus.GetA_asInt() >= 1 &&
             aTrack.GetKineticEnergy() > 2340*CLHEP::MeV )    
  {
    applied = true;
  }
  else if( ( aTrack.GetDefinition() == G4KaonPlus::KaonPlus() || 
	     aTrack.GetDefinition() == G4KaonMinus::KaonMinus()  ) && 
             targetNucleus.GetA_asInt() >= 1 &&
             aTrack.GetKineticEnergy() > 1980*CLHEP::MeV ) 
  {
    applied = true;
  }
  return applied;
}

/////////////////////////////////////////////////////////////////
//
// Return dissociated particle products and recoil nucleus

G4HadFinalState*  
G4LMsdGenerator::ApplyYourself( const G4HadProjectile& aTrack, 
                                      G4Nucleus& targetNucleus )
{
  theParticleChange.Clear();

  const G4HadProjectile* aParticle = &aTrack;
  G4double eTkin = aParticle->GetKineticEnergy();

  if( eTkin <= 1.*CLHEP::GeV && aTrack.GetDefinition() != G4Proton::Proton()) 
  {
    theParticleChange.SetEnergyChange(eTkin);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

  G4int A = targetNucleus.GetA_asInt();
  G4int Z = targetNucleus.GetZ_asInt();

  G4double plab = aParticle->GetTotalMomentum();
  G4double plab2 = plab*plab;
  
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double partMass = theParticle->GetPDGMass();

  G4double oldE = partMass + eTkin;

  G4double targMass   = G4NucleiProperties::GetNuclearMass(A, Z);
  G4double targMass2   = targMass*targMass;

  G4LorentzVector partLV = aParticle->Get4Momentum();

  G4double sumE = oldE + targMass;
  G4double sumE2 = sumE*sumE;
   
  G4ThreeVector p1     = partLV.vect();
  // G4cout<<"p1 = "<<p1<<G4endl;
  G4ParticleMomentum p1unit = p1.unit();

  G4double Mx = SampleMx(aParticle); // in GeV
  G4double t  = SampleT( aParticle, Mx); // in GeV

  Mx *= CLHEP::GeV;

  G4double Mx2 = Mx*Mx;

  // equation for q|| based on sum-E-P and new invariant mass

  G4double B = sumE2 + targMass2 - Mx2 - plab2;

  G4double a = 4*(plab2 - sumE2);
  G4double b = 4*plab*B;
  G4double c = B*B - 4*sumE2*targMass2;
  G4double det2 = b*b - 4*a*c;
  G4double qLong,  det, eRetard; // , x2, x3, e2;

  if( det2 >= 0.) 
  {
    det = std::sqrt(det2);
    qLong = (-b - det)/2./a;
    eRetard = std::sqrt((plab-qLong)*(plab-qLong)+Mx2);
  }
  else 
  {
    theParticleChange.SetEnergyChange(eTkin);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }
  theParticleChange.SetStatusChange(stopAndKill);

  plab -= qLong;

  G4ThreeVector pRetard = plab*p1unit;

  G4ThreeVector pTarg = p1 - pRetard; 

  G4double eTarg = std::sqrt( targMass2 + pTarg.mag2()); // std::sqrt( targMass*targMass + pTarg.mag2() );

  G4LorentzVector lvRetard(pRetard, eRetard);
  G4LorentzVector lvTarg(pTarg, eTarg);
   
  lvTarg += lvRetard; // sum LV

  G4ThreeVector bst = lvTarg.boostVector();

  lvRetard.boost(-bst); // to CNS

  G4ThreeVector pCMS   = lvRetard.vect();
  G4double momentumCMS = pCMS.mag();
  G4double tMax        = 4.0*momentumCMS*momentumCMS;

  if( t > tMax ) t = tMax*G4UniformRand();

  G4double cost = 1. - 2.0*t/tMax;
  

  G4double phi  = G4UniformRand()*CLHEP::twopi;
  G4double sint;

  if( cost > 1.0 || cost < -1.0 ) // 
  {
    cost = 1.0;
    sint = 0.0;    
  } 
  else  // normal situation
  {
    sint = std::sqrt( (1.0-cost)*(1.0+cost) );
  }    
  G4ThreeVector v1( sint*std::cos(phi), sint*std::sin(phi), cost);

  v1 *= momentumCMS;
  
  G4LorentzVector lvRes( v1.x(),v1.y(),v1.z(), std::sqrt( momentumCMS*momentumCMS + Mx2));

  lvRes.boost(bst); // to LS

  lvTarg -= lvRes;

  G4double eRecoil =  lvTarg.e() - targMass;

  if( eRecoil > 100.*CLHEP::MeV ) // add recoil nucleus
  {
    G4ParticleDefinition * recoilDef = 0;

    if      ( Z == 1 && A == 1 ) { recoilDef = G4Proton::Proton(); }
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

  G4ParticleDefinition* ddPart = G4ParticleTable::GetParticleTable()->
                                 FindParticle(fPDGencoding);

  // G4cout<<fPDGencoding<<", "<<ddPart->GetParticleName()<<", "<<ddPart->GetPDGMass()<<" MeV; lvRes = "<<lvRes<<G4endl;

  // G4DynamicParticle * aRes = new G4DynamicParticle( ddPart, lvRes);
  //  theParticleChange.AddSecondary(aRes); // simply return resonance


  
  // Recursive decay using methods of G4KineticTrack

  G4KineticTrack ddkt( ddPart, 0., G4ThreeVector(0.,0.,0.), lvRes);
  G4KineticTrackVector* ddktv = ddkt.Decay();

  G4DecayKineticTracks decay( ddktv );

  for( unsigned int i = 0; i < ddktv->size(); i++ ) // add products to partchange
  {
    G4DynamicParticle * aNew = 
      new G4DynamicParticle( ddktv->operator[](i)->GetDefinition(),
                             ddktv->operator[](i)->Get4Momentum());

    // G4cout<<"       "<<i<<", "<<aNew->GetDefinition()->GetParticleName()<<", "<<aNew->Get4Momentum()<<G4endl;

    theParticleChange.AddSecondary(aNew);
    delete ddktv->operator[](i);
  }
  delete ddktv;
   
  return &theParticleChange;
}

//////////////////////////////////////
//
// Sample Mx as Roper resonances, set PDG encoding 

G4double G4LMsdGenerator::SampleMx(const G4HadProjectile* aParticle)                          
{
  G4double Mx = 0.;
  G4int i;
  G4double rand = G4UniformRand();

  for( i = 0; i < 60; i++)
  {
    if( rand >= fProbMx[i][1] ) break;
  }
  if(i <= 0)       Mx = fProbMx[0][0];
  else if(i >= 59) Mx = fProbMx[59][0];
  else             Mx = fProbMx[i][0];

  fPDGencoding = 0;

  if ( Mx <=  1.45 )   
  {   
    if( aParticle->GetDefinition() == G4Proton::Proton() )
    {        
      Mx = 1.44;
      // fPDGencoding = 12212;
      fPDGencoding = 2214;
    }
    else if( aParticle->GetDefinition() == G4Neutron::Neutron() ) 
    {
      Mx = 1.44;
      fPDGencoding = 12112;
    }
    else if( aParticle->GetDefinition() == G4PionPlus::PionPlus() ) 
    {
      // Mx = 1.3;
      // fPDGencoding = 100211;
      Mx = 1.26;
      fPDGencoding = 20213; // a1(1260)+
    }
    else if( aParticle->GetDefinition() == G4PionMinus::PionMinus() ) 
    {
      // Mx = 1.3;
      // fPDGencoding = -100211;
      Mx = 1.26;
      fPDGencoding = -20213; // a1(1260)-
    }
    else if( aParticle->GetDefinition() == G4KaonPlus::KaonPlus() ) 
    {
      Mx = 1.27;
      fPDGencoding = 10323;
    }
    else if( aParticle->GetDefinition() == G4KaonMinus::KaonMinus() ) 
    {
      Mx = 1.27;
      fPDGencoding = -10323;
    }
  }
  else if ( Mx <=  1.55 ) 
  {   
    if( aParticle->GetDefinition() == G4Proton::Proton() ) 
    {       
      Mx = 1.52;
      // fPDGencoding = 2124;
      fPDGencoding = 2214;
    }
    else if( aParticle->GetDefinition() == G4Neutron::Neutron() ) 
    {
      Mx = 1.52;
      fPDGencoding = 1214;
    }
    else if( aParticle->GetDefinition() == G4PionPlus::PionPlus() ) 
    {
      // Mx = 1.45;
      // fPDGencoding = 10211;
      Mx = 1.32;
      fPDGencoding = 215; // a2(1320)+
    }
    else if( aParticle->GetDefinition() == G4PionMinus::PionMinus() ) 
    {
      // Mx = 1.45;
      // fPDGencoding = -10211;
      Mx = 1.32;
      fPDGencoding = -215; // a2(1320)-
    }
    else if( aParticle->GetDefinition() == G4KaonPlus::KaonPlus() ) 
    {
      Mx = 1.46;
      fPDGencoding = 100321;
    }
    else if( aParticle->GetDefinition() == G4KaonMinus::KaonMinus() ) 
    {
      Mx = 1.46;
      fPDGencoding = -100321;
    }
  }
  else                    
  {    
    if( aParticle->GetDefinition() == G4Proton::Proton() ) 
    {       
      Mx = 1.68;
      // fPDGencoding = 12216;
      fPDGencoding = 2214;
    }
    else if( aParticle->GetDefinition() == G4Neutron::Neutron() ) 
    {
      Mx = 1.68;
      fPDGencoding = 12116;
    }
    else if( aParticle->GetDefinition() == G4PionPlus::PionPlus() ) 
    {
      Mx = 1.67;
      fPDGencoding = 10215;  // pi2(1670)+
      // Mx = 1.45;
      // fPDGencoding = 10211;
    }
    else if( aParticle->GetDefinition() == G4PionMinus::PionMinus() ) 
    {
      Mx = 1.67;     // f0 problems->4pi vmg 20.11.14
      fPDGencoding = -10215;  // pi2(1670)-
      // Mx = 1.45;
      // fPDGencoding = -10211;
    }
    else if( aParticle->GetDefinition() == G4KaonPlus::KaonPlus() ) 
    {
      Mx = 1.68;
      fPDGencoding = 30323;
    }
    else if( aParticle->GetDefinition() == G4KaonMinus::KaonMinus() ) 
    {
      Mx = 1.68;
      fPDGencoding = -30323;
    }
  } 
  if(fPDGencoding == 0)
  {
      Mx = 1.44;
      // fPDGencoding = 12212;   
      fPDGencoding = 2214;
  } 
  G4ParticleDefinition* myResonance = 
  G4ParticleTable::GetParticleTable()->FindParticle( fPDGencoding );
  
  if ( myResonance ) Mx = myResonance->GetPDGMass();

  // G4cout<<"PDG-ID = "<<fPDGencoding<<"; with mass = "<<Mx/CLHEP::GeV<<" GeV"<<G4endl;

  return Mx/CLHEP::GeV;
}

////////////////////////////////////// 
//
// Sample t with kinematic limitations of Mx and Tkin

G4double G4LMsdGenerator::SampleT(  const G4HadProjectile* aParticle, 
G4double Mx)                          
{
  G4double t=0., b=0., rTkin = 50.*CLHEP::GeV, eTkin = aParticle->GetKineticEnergy();
  G4int i;

  for( i = 0; i < 23; ++i)
  {
    if( Mx <= fMxBdata[i][0] ) break;
  }
  if( i <= 0 )       b = fMxBdata[0][1];
  else if( i >= 22 ) b = fMxBdata[22][1];
  else               b = fMxBdata[i][1];

  if( eTkin > rTkin ) b *= 1. + G4Log(eTkin/rTkin);

  G4double rand = G4UniformRand();

  t = -G4Log(rand)/b;

  t *= (CLHEP::GeV*CLHEP::GeV); // in G4 internal units

  return t;
}


////////////////////////////////////////////////
//
// Integral spectrum of Mx (GeV)

const G4double G4LMsdGenerator::fProbMx[60][2] = 
{
  {1.000000e+00, 	1.000000e+00},
  {1.025000e+00, 	1.000000e+00},
  {1.050000e+00, 	1.000000e+00},
  {1.075000e+00, 	1.000000e+00},
  {1.100000e+00, 	9.975067e-01},
  {1.125000e+00, 	9.934020e-01},
  {1.150000e+00, 	9.878333e-01},
  {1.175000e+00, 	9.805002e-01},
  {1.200000e+00, 	9.716846e-01},
  {1.225000e+00, 	9.604761e-01},
  {1.250000e+00, 	9.452960e-01},
  {1.275000e+00, 	9.265278e-01},
  {1.300000e+00, 	9.053632e-01},
  {1.325000e+00, 	8.775566e-01},
  {1.350000e+00, 	8.441969e-01},
  {1.375000e+00, 	8.076336e-01},
  {1.400000e+00, 	7.682520e-01},
  {1.425000e+00, 	7.238306e-01},
  {1.450000e+00, 	6.769306e-01},
  {1.475000e+00, 	6.303898e-01},
  {1.500000e+00, 	5.824632e-01},
  {1.525000e+00, 	5.340696e-01},
  {1.550000e+00, 	4.873736e-01},
  {1.575000e+00, 	4.422901e-01},
  {1.600000e+00, 	3.988443e-01},
  {1.625000e+00, 	3.583727e-01},
  {1.650000e+00, 	3.205405e-01},
  {1.675000e+00, 	2.856655e-01},
  {1.700000e+00, 	2.537508e-01},
  {1.725000e+00, 	2.247863e-01},
  {1.750000e+00, 	1.985798e-01},
  {1.775000e+00, 	1.750252e-01},
  {1.800000e+00, 	1.539777e-01},
  {1.825000e+00, 	1.352741e-01},
  {1.850000e+00, 	1.187157e-01},
  {1.875000e+00, 	1.040918e-01},
  {1.900000e+00, 	9.118422e-02},
  {1.925000e+00, 	7.980909e-02},
  {1.950000e+00, 	6.979378e-02},
  {1.975000e+00, 	6.097771e-02},
  {2.000000e+00, 	5.322122e-02},
  {2.025000e+00, 	4.639628e-02},
  {2.050000e+00, 	4.039012e-02},
  {2.075000e+00, 	3.510275e-02},
  {2.100000e+00, 	3.044533e-02},
  {2.125000e+00, 	2.633929e-02},
  {2.150000e+00, 	2.271542e-02},
  {2.175000e+00, 	1.951295e-02},
  {2.200000e+00, 	1.667873e-02},
  {2.225000e+00, 	1.416633e-02},
  {2.250000e+00, 	1.193533e-02},
  {2.275000e+00, 	9.950570e-03},
  {2.300000e+00, 	8.181515e-03},
  {2.325000e+00, 	6.601664e-03},
  {2.350000e+00, 	5.188025e-03},
  {2.375000e+00, 	3.920655e-03},
  {2.400000e+00, 	2.782246e-03},
  {2.425000e+00, 	1.757765e-03},
  {2.450000e+00, 	8.341435e-04},
  {2.475000e+00, 	0.000000e+00}
};

//////////////////////////////////////////////
//
// Slope b (1/GeV/GeV) vs Mx (GeV) for t-sampling over exp(-b*t) 

const G4double G4LMsdGenerator::fMxBdata[23][2] = 
{
  {1.09014,      17.8620},      
  {1.12590,      19.2831},    
  {1.18549,      17.6907},      
  {1.21693,      16.4760},      
  {1.25194,      15.3867},      
  {1.26932,      14.4236},      
  {1.29019,      13.2931},      
  {1.30755,      12.2882},      
  {1.31790,      11.4509},      
  {1.33888,      10.6969},      
  {1.34911,      9.44130},      
  {1.37711,      8.56148},      
  {1.39101,      7.76593},      
  {1.42608,      6.88582},      
  {1.48593,      6.13019},      
  {1.53179,      5.87723},      
  {1.58111,      5.37308},      
  {1.64105,      4.95217},      
  {1.69037,      4.44803},      
  {1.81742,      3.89879},      
  {1.88096,      3.68693},      
  {1.95509,      3.43278},      
  {2.02219,      3.30445}
};



//
//
/////////////////////////////////////////////
