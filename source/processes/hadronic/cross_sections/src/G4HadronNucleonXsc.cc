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
// 14.03.07 V. Grichine - first implementation
//
// 04.09.18 V. Ivantchenko Major revision of interfaces and implementation
// 30.09.18 V. Grichine hyperon-nucleon xsc first implementation

#include "G4HadronNucleonXsc.hh"
#include "G4Element.hh"
#include "G4Proton.hh"
#include "G4Nucleus.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"
#include "G4NuclearRadii.hh"

#include "G4LambdacPlus.hh"
#include "G4AntiLambdacPlus.hh"
#include "G4AntiXibZero.hh"
#include "G4OmegacZero.hh"
#include "G4SigmacZero.hh"
#include "G4AntiLambdab.hh"
#include "G4AntiSigmabMinus.hh"
#include "G4AntiXicPlus.hh"
#include "G4AntiLambdacPlus.hh"
#include "G4AntiSigmabPlus.hh"
#include "G4AntiXicZero.hh"
#include "G4AntiSigmabZero.hh"
#include "G4XibMinus.hh"
#include "G4AntiSigmacPlus.hh"
#include "G4XibZero.hh"
#include "G4AntiOmegabMinus.hh"
#include "G4AntiSigmacPlusPlus.hh"

#include "G4Lambdab.hh"
#include "G4SigmabMinus.hh"
#include "G4XicPlus.hh"
#include "G4AntiOmegacZero.hh"
#include "G4AntiSigmacZero.hh"
#include "G4LambdacPlus.hh"
#include "G4SigmabPlus.hh"
#include "G4XicZero.hh"
#include "G4SigmabZero.hh"
#include "G4SigmacPlus.hh"
#include "G4AntiXibMinus.hh"
#include "G4OmegabMinus.hh"
#include "G4SigmacPlusPlus.hh"

#include "G4BMesonZero.hh"
#include "G4AntiBMesonZero.hh"
#include "G4DMesonZero.hh"
#include "G4AntiDMesonZero.hh"
#include "G4BsMesonZero.hh"
#include "G4AntiBsMesonZero.hh"
#include "G4BcMesonPlus.hh"
#include "G4BcMesonMinus.hh"
#include "G4DsMesonPlus.hh"
#include "G4DsMesonMinus.hh"
#include "G4Eta.hh"
#include "G4EtaPrime.hh"
#include "G4Etac.hh"

#include "G4BMesonPlus.hh"
#include "G4BMesonMinus.hh"

#include "G4DMesonPlus.hh"
#include "G4DMesonMinus.hh"

#include "G4JPsi.hh"
#include "G4Upsilon.hh"

#include "G4Lambda.hh"
#include "G4AntiLambda.hh"
#include "G4SigmaPlus.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4SigmaMinus.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4SigmaZero.hh"
#include "G4AntiSigmaZero.hh"
#include "G4XiMinus.hh"
#include "G4XiZero.hh"
#include "G4AntiXiMinus.hh"
#include "G4AntiXiZero.hh"
#include "G4OmegaMinus.hh"
#include "G4AntiOmegaMinus.hh"

static const G4double invGeV  = 1.0/CLHEP::GeV;
static const G4double invGeV2 = 1.0/(CLHEP::GeV*CLHEP::GeV);
// PDG fit constants
static const G4double minLogP = 3.5;    // min of (lnP-minLogP)^2 
static const G4double cofLogE = .0557;  // elastic (lnP-minLogP)^2 
static const G4double cofLogT = .3;     // total (lnP-minLogP)^2 
static const G4double pMin = .1;        // fast LE calculation 
static const G4double pMax = 1000.;     // fast HE calculation 
static const G4double ekinmin = 0.1*CLHEP::MeV; // protection against zero ekin 

G4HadronNucleonXsc::G4HadronNucleonXsc() 
  : fTotalXsc(0.0), fElasticXsc(0.0), fInelasticXsc(0.0)
{
  fHypTotXscCof = 0.88; // for transformation pp(pn) to hyperon-nucleon

  theGamma    = G4Gamma::Gamma();
  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theAProton  = G4AntiProton::AntiProton();
  theANeutron = G4AntiNeutron::AntiNeutron();
  thePiPlus   = G4PionPlus::PionPlus();
  thePiMinus  = G4PionMinus::PionMinus();
  thePiZero   = G4PionZero::PionZero();
  theD        = G4Deuteron::Deuteron();
  theT        = G4Triton::Triton();
  theA        = G4Alpha::Alpha();
  theHe3      = G4He3::He3();
  // strange
  theKPlus    = G4KaonPlus::KaonPlus();
  theKMinus   = G4KaonMinus::KaonMinus();
  theK0S      = G4KaonZeroShort::KaonZeroShort();
  theK0L      = G4KaonZeroLong::KaonZeroLong();
  theL        = G4Lambda::Lambda();
  theAntiL    = G4AntiLambda::AntiLambda();
  theSPlus    = G4SigmaPlus::SigmaPlus();
  theASPlus   = G4AntiSigmaPlus::AntiSigmaPlus();
  theSMinus   = G4SigmaMinus::SigmaMinus();
  theASMinus  = G4AntiSigmaMinus::AntiSigmaMinus();
  theS0       = G4SigmaZero::SigmaZero();
  theAS0      = G4AntiSigmaZero::AntiSigmaZero();
  theXiMinus  = G4XiMinus::XiMinus();
  theXi0      = G4XiZero::XiZero();
  theAXiMinus = G4AntiXiMinus::AntiXiMinus();
  theAXi0     = G4AntiXiZero::AntiXiZero();
  theOmega    = G4OmegaMinus::OmegaMinus();
  theAOmega   = G4AntiOmegaMinus::AntiOmegaMinus();
  // c- and b- hyperons 
  theLambdaCPlus = G4LambdacPlus::LambdacPlus();
  theALambdaCPlus = G4AntiLambdacPlus::AntiLambdacPlus();
  theOmegaC0 = G4OmegacZero::OmegacZero();
  theAOmegaC0 = G4AntiOmegacZero::AntiOmegacZero();
  theSigmaCPlus = G4SigmacPlus::SigmacPlus();
  theASigmaCPlus = G4AntiSigmacPlus::AntiSigmacPlus();
  theSigmacPP =  G4SigmacPlusPlus::SigmacPlusPlus();
  theASigmacPP =  G4AntiSigmacPlusPlus::AntiSigmacPlusPlus();
  theSigmaC0 = G4SigmacZero::SigmacZero();
  theASigmaC0 = G4AntiSigmacZero::AntiSigmacZero();
  theXiCPlus = G4XicPlus::XicPlus();
  theAXiCPlus = G4AntiXicPlus::AntiXicPlus();
  theXiC0 = G4XicZero::XicZero();
  theAXiC0 = G4AntiXicZero::AntiXicZero();
  theLambdaB = G4Lambdab::Lambdab();
  theALambdaB = G4AntiLambdab::AntiLambdab();
  theOmegaBMinus = G4OmegabMinus::OmegabMinus();
  theAOmegaBMinus = G4AntiOmegabMinus::AntiOmegabMinus();
  theSigmaBMinus = G4SigmabMinus::SigmabMinus();
  theASigmaBMinus = G4AntiSigmabMinus::AntiSigmabMinus();
  theSigmaBPlus = G4SigmabPlus::SigmabPlus();
  theASigmaBPlus = G4AntiSigmabPlus::AntiSigmabPlus();
  theSigmaB0 = G4SigmabZero::SigmabZero();
  theASigmaB0 = G4AntiSigmabZero::AntiSigmabZero();
  theXiBMinus = G4XibMinus::XibMinus();
  theAXiBMinus = G4AntiXibMinus::AntiXibMinus();
  theXiB0 = G4XibZero::XibZero();
  theAXiB0 = G4AntiXibZero::AntiXibZero();
  //(s-) c- and b-mesons
  theBMeson0 = G4BMesonZero::BMesonZero();
  theABMeson0 = G4AntiBMesonZero::AntiBMesonZero();
  theDMeson0 = G4DMesonZero::DMesonZero();
  theADMeson0 = G4AntiDMesonZero::AntiDMesonZero();
  theBsMeson0 = G4BsMesonZero::BsMesonZero();
  theABsMeson0 = G4AntiBsMesonZero::AntiBsMesonZero();
  theBcMesonPlus = G4BcMesonPlus::BcMesonPlus();
  theBcMesonMinus = G4BcMesonMinus::BcMesonMinus();
  theDsMesonPlus = G4DsMesonPlus::DsMesonPlus();
  theDsMesonMinus = G4DsMesonMinus::DsMesonMinus();
  theDMesonPlus = G4DMesonPlus::DMesonPlus();
  theDMesonMinus = G4DMesonMinus::DMesonMinus();
  theBMesonPlus = G4BMesonPlus::BMesonPlus();
  theBMesonMinus = G4BMesonMinus::BMesonMinus();
  theEta = G4Eta::Eta();
  theEtaPrime = G4EtaPrime::EtaPrime();
  theEtaC = G4Etac::Etac();
  theJPsi = G4JPsi::JPsi();
  theUpsilon = G4Upsilon::Upsilon();

  g4calc = G4Pow::GetInstance();
}

G4HadronNucleonXsc::~G4HadronNucleonXsc()
{}

void G4HadronNucleonXsc::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4HadronNucleonXsc calculates the total, inelastic and elastic\n"
          << "hadron-nucleon cross sections using several different\n"
          << "parameterizations within the Glauber-Gribov approach. It is\n"
          << "valid for all incident gammas and long-lived hadrons at\n"
          << "energies above 30 keV.  This is a cross section component which\n"
          << "is to be used to build a cross section data set.\n"; 
}



G4double G4HadronNucleonXsc::HadronNucleonXsc(  const G4ParticleDefinition* theParticle,
                               const G4ParticleDefinition* nucleon, G4double ekin)
{
  G4double xsc(0.);
  G4int    pdg = std::abs( theParticle->GetPDGEncoding() );
  
  if ( pdg == 2212 || pdg == 2112 || pdg == 211 ) // p, n, pi+-
  {
    xsc = HadronNucleonXscNS( theParticle, nucleon, ekin);
  }
  else if ( pdg == 321 || pdg == 310 || pdg == 130 ) // K+-, K0, Ks
  {
    xsc = KaonNucleonXscNS( theParticle, nucleon, ekin);
  }
  else if ( pdg == 3122 || pdg == 3222 || pdg == 3112 || pdg == 3212 || pdg == 3322 || pdg == 3312 || pdg == 3324 ||

	    pdg == 4122 || pdg == 4332 || pdg == 4122 || pdg == 4212 || pdg == 4222 || pdg == 4112 || pdg == 4232 || pdg == 4132  || 

            pdg == 5122 || pdg == 5332 || pdg == 5122 || pdg == 5112 || pdg == 5222 || pdg == 5212 || pdg == 5132 || pdg == 5232  
	    
	  ) // heavy s-,c-,b-hyperons
  {
    xsc = HyperonNucleonXscNS( theParticle, nucleon, ekin);
  }
  else if ( pdg == 511 || pdg == 421 || pdg == 531 || pdg == 541 || pdg == 431 || pdg == 411 || pdg == 521 ||

	    pdg == 221 || pdg == 331 || pdg == 441 || pdg == 443 || pdg == 543 

	    ) // s-,c-,b-mesons
  {
    xsc = SCBMesonNucleonXscNS( theParticle, nucleon, ekin);
  }
  else
  {
    xsc = HadronNucleonXscNS( theParticle, nucleon, ekin);
  }
  return xsc;
}


//////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to PDG parametrisation (2017):
// http://pdg.lbl.gov/2017/reviews/hadronicrpp.pdf

G4double G4HadronNucleonXsc::HadronNucleonXscPDG(
                             const G4ParticleDefinition* theParticle, 
			     const G4ParticleDefinition* nucleon, G4double ekin)
{
  static const G4double M    = 2.1206; // in Gev
  static const G4double eta1 = 0.4473;
  static const G4double eta2 = 0.5486;
  static const G4double H    = 0.272;

  G4double mass1 = theParticle->GetPDGMass();
  if(theParticle == theGamma) { mass1 = 770.; }
  G4double mass2 = nucleon->GetPDGMass();

  G4double sMand = CalcMandelstamS(ekin, mass1, mass2)*invGeV2;
  G4double x = (mass1 + mass2)*invGeV + M;
  G4double blog = G4Log(sMand/(x*x));

  G4double P(0.0), R1(0.0), R2(0.0), del(1.0);

  G4bool proton  = (nucleon == theProton);
  G4bool neutron = (nucleon == theNeutron);
  
  if(theParticle == theNeutron)
  {
    if ( proton )
    {
      P  = 34.71;
      R1 = 12.52;
      R2 = -6.66;
    }
    else
    {
      P  = 34.41;
      R1 = 13.07;
      R2 = -7.394;
    }
  } 
  else if(theParticle == theProton) 
  {
    if ( neutron )
    {
      P  = 34.71;
      R1 = 12.52;
      R2 = -6.66;
    }
    else
    {
      P  = 34.41;
      R1 = 13.07;
      R2 = -7.394;
    }
  } 
  else if(theParticle == theAProton) 
  {
    if ( neutron )
    {
      P  = 34.71;
      R1 = 12.52;
      R2 = 6.66;
    }
    else
    {
      P  = 34.41;
      R1 = 13.07;
      R2 = 7.394;
    }
  } 
  else if(theParticle == theANeutron) 
  {
    if ( proton )
    {
      P  = 34.71;
      R1 = 12.52;
      R2 = 6.66;
    }
    else
    {
      P  = 34.41;
      R1 = 13.07;
      R2 = 7.394;
    }
  } 
  else if(theParticle == thePiPlus) 
  {
    P  = 18.75;
    R1 = 9.56;
    R2 = -1.767;
  } 
  else if(theParticle == thePiMinus) 
  {
    P  = 18.75;
    R1 = 9.56;
    R2 = 1.767;
  } 
  else if(theParticle == theKPlus) 
  {
    if ( proton )
    {
      P  = 16.36;
      R1 = 4.29;
      R2 = -3.408;
    }
    else
    {
      P  = 16.31;
      R1 = 3.7;
      R2 = -1.826;
    }
  } 
  else if(theParticle == theKMinus) 
  {
    if ( proton )
    {
      P  = 16.36;
      R1 = 4.29;
      R2 = 3.408;
    }
    else
    {
      P  = 16.31;
      R1 = 3.7;
      R2 = 1.826;
    }
  }
  else if(theParticle == theK0S || theParticle == theK0L) 
  {
    P  = 16.36;
    R1 = 2.5;
    R2 = 0.;
  }
  else if(theParticle == theSMinus) 
  {
    P  = 34.7;
    R1 = -46.;
    R2 = 48.;
  } 
  else if(theParticle == theGamma) // modify later on
  {
    del= 0.003063;
    P  = 34.71*del;
    R1 = (neutron) ? 0.0231 : 0.0139;
    R2 = 0.;
     
  } 
  else  // as proton ??? 
  {
    if ( neutron )
    {
      P  = 34.71;
      R1 = 12.52;
      R2 = -6.66;
    }
    else
    {
      P  = 34.41;
      R1 = 13.07;
      R2 = -7.394;
    }
  } 
  fTotalXsc = CLHEP::millibarn* 
    (del*(H*blog*blog + P) + R1*G4Exp(-eta1*blog) + R2*G4Exp(-eta2*blog));
  fInelasticXsc = 0.75*fTotalXsc;
  fElasticXsc   = fTotalXsc - fInelasticXsc;

  if( proton && theParticle->GetPDGCharge() > 0. && ekin < 100*MeV)
  {
    G4double cB = CoulombBarrier(theParticle, nucleon, ekin);
    fTotalXsc   *= cB;
    fElasticXsc *= cB; 
    fInelasticXsc *= cB; 
  }
  /*
  G4cout << "## HadronNucleonXscPDG: ekin(MeV)= " << ekin 
	 << " tot= " << fTotalXsc/CLHEP::millibarn 
	 << " inel= " << fInelasticXsc/CLHEP::millibarn
	 << " el= " << fElasticXsc/CLHEP::millibarn
         << G4endl;
  */
  return fTotalXsc;
}

//////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon cross-section based on N. Starkov parametrisation of
// data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database

G4double G4HadronNucleonXsc::HadronNucleonXscNS(
                             const G4ParticleDefinition* theParticle, 
			     const G4ParticleDefinition* nucleon, G4double ekin0)
{  
  const G4double ekin = std::max(ekin0, ekinmin);
  /*  
  G4cout<< "HadronNucleonXscNS: Ekin(GeV)= " << ekin/GeV << "  "
        << theParticle->GetParticleName() << " + " 
        << nucleon->GetParticleName() << G4endl;
  */
  if(theParticle == theAProton || theParticle == theANeutron) {
    return HadronNucleonXscPDG(theParticle, nucleon, ekin);
  }

  G4double pM   = theParticle->GetPDGMass();
  G4double tM   = nucleon->GetPDGMass();
  G4double pE   = ekin + pM; 
  G4double pLab = std::sqrt(ekin*(ekin + 2*pM));

  G4double sMand = CalcMandelstamS(ekin, pM, tM)*invGeV2;

  pLab *= invGeV;
  pE   *= invGeV;

  if(pLab >= 10.) {
    fTotalXsc = HadronNucleonXscPDG(theParticle, nucleon, ekin)/CLHEP::millibarn;
  } else { fTotalXsc = 0.0; }
  fElasticXsc = 0.0;
  //G4cout << "Stot(mb)= " << fTotalXsc << " pLab= " << pLab 
  //	 << " Smand= " << sMand <<G4endl;
  G4double logP = G4Log(pLab);

  G4bool proton = (nucleon == theProton);
  G4bool neutron = (nucleon == theNeutron);

  if( theParticle == theNeutron) 
  {
    if( pLab >= 373.)
    {
      fElasticXsc = 6.5 + 0.308*G4Exp(G4Log(G4Log(sMand/400.)*1.65)) 
	+ 9.19*G4Exp(-G4Log(sMand)*0.458);
    }
    else if( pLab >= 100.)
    {
      fElasticXsc = 5.53 + 0.308*G4Exp(G4Log(G4Log(sMand/28.9))*1.1) 
	+ 9.19*G4Exp(-G4Log(sMand)*0.458);
    }
    else if( pLab >= 10.)
    {
      fElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
    }
    else  // pLab < 10 GeV/c
    {
      if( neutron )      // nn to be pp
      {
        G4double x = G4Log(pLab/0.73);
        if( pLab < 0.4 )
        {
          fTotalXsc = 23 + 50*std::sqrt(g4calc->powN(-x, 7));
          fElasticXsc = fTotalXsc;
        }
        else if( pLab < 0.73 )
        {
          fTotalXsc = 23 + 50*std::sqrt(g4calc->powN(-x, 7));
          fElasticXsc = fTotalXsc; 
        }
        else if( pLab < 1.05  )
        {
          fTotalXsc = 23 + 40*x*x;
          fElasticXsc = 23 + 20*x*x;
        }
        else    // 1.05 - 10 GeV/c
        {
          fTotalXsc = 39.0+75*(pLab - 1.2)/(g4calc->powN(pLab,3) + 0.15);
          fElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
        }
      }
      if( proton )   // pn to be np
      {
        if( pLab < 0.02 )
        {
          fTotalXsc = 4100+30*G4Exp(G4Log(G4Log(1.3/pLab))*3.6); 
	  fElasticXsc = fTotalXsc;
        }      
        else if( pLab < 0.8 )
        {
          fTotalXsc = 33+30*g4calc->powN(G4Log(pLab/1.3),4);
	  fElasticXsc = fTotalXsc;
        }      
        else if( pLab < 1.4 )
        {
          fTotalXsc = 33+30*g4calc->powN(G4Log(pLab/0.95),2);
          G4double x = G4Log(0.511/pLab);
          fElasticXsc =  6 + 52/( x*x + 1.6 );
        }
        else    // 1.4 < pLab < 10.  )
        {
          fTotalXsc = 33.3 + 20.8*(pLab*pLab - 1.35)
	    /(std::sqrt(g4calc->powN(pLab,5)) + 0.95);
          fElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
        }
      }
    }
  }
  ////// proton //////////////////////////////////////////////
  else if(theParticle == theProton) 
  {
    if( pLab >= 373.) // pdg due to TOTEM data
    {
      fElasticXsc = 6.5 + 0.308*G4Exp(G4Log(G4Log(sMand/400.))*1.65) 
	+ 9.19*G4Exp(-G4Log(sMand)*0.458);     
    }
    else if( pLab >= 100.)
    {
      fElasticXsc = 5.53 + 0.308*G4Exp(G4Log(G4Log(sMand/28.9))*1.1) 
	+ 9.19*G4Exp(-G4Log(sMand)*0.458);     
    }
    else if( pLab >= 10.)
    {
      fElasticXsc = 6. + 20./( (logP-0.182)*(logP-0.182) + 1.0 );
    }
    else
    {
      // pp
      if( proton )
      {
        if( pLab < 0.73 )
        {
          fTotalXsc = 23 + 50*std::sqrt(g4calc->powN(G4Log(0.73/pLab),7));
          fElasticXsc = fTotalXsc;
        }
        else if( pLab < 1.05  )
        {
	  G4double x = G4Log(pLab/0.73);
          fTotalXsc = 23 + 40*x*x;
          fElasticXsc = 23 + 20*x*x;
        }
        else    // 1.05 - 10 GeV/c
        {
          fTotalXsc = 39.0+75*(pLab - 1.2)/(g4calc->powN(pLab,3) + 0.15);
          fElasticXsc = 6. + 20./( (logP-0.182)*(logP-0.182) + 1.0 );
        }
      }
      else if( neutron )     // pn to be np
      {
        if( pLab < 0.02 )
        {
          fTotalXsc = 4100+30*G4Exp(G4Log(G4Log(1.3/pLab))*3.6); 
	  fElasticXsc = fTotalXsc;
        }      
        else if( pLab < 0.8 )
        {
          fTotalXsc = 33+30*g4calc->powN(G4Log(pLab/1.3),4);
	  fElasticXsc = fTotalXsc;
        }      
        else if( pLab < 1.4 )
        {
          G4double x1 = G4Log(pLab/0.95);
          G4double x2 = G4Log(0.511/pLab);
          fTotalXsc = 33+30*x1*x1;
          fElasticXsc =  6 + 52/(x2*x2 + 1.6);
        }
        else    // 1.4 < pLab < 10.  )
        {
          fTotalXsc = 33.3 + 20.8*(pLab*pLab - 1.35)
	    /(std::sqrt(g4calc->powN(pLab,5)) + 0.95);
          fElasticXsc = 6. + 20./( (logP-0.182)*(logP-0.182) + 1.0 );
        }
      }
    }    
  }
  // pi+ p; pi- n 
  else if((theParticle == thePiPlus && proton) || 
	  (theParticle == thePiMinus && neutron))  
  {
    if( pLab < 0.28 )
    {
      fTotalXsc = 10./((logP + 1.273)*(logP + 1.273) + 0.05);
      fElasticXsc = fTotalXsc;
    }
    else if( pLab < 0.68 )
    {
      fTotalXsc = 14./((logP + 1.273)*(logP + 1.273) + 0.07);
      fElasticXsc = fTotalXsc;
    }
    else if( pLab < 0.85 )
    {
      G4double x = G4Log(pLab/0.77);
      fTotalXsc = 88.*x*x + 14.9;
      fElasticXsc = fTotalXsc*G4Exp(-3.*(pLab - 0.68));  
    }
    else if( pLab < 1.15 )
    {
      G4double x = G4Log(pLab/0.77);
      fTotalXsc = 88.*x*x + 14.9;
      fElasticXsc = 6.0 + 1.4/(( pLab - 1.4)*( pLab - 1.4) + 0.1);
    }
    else if( pLab < 1.4) // ns original
    {
      G4double Ex1 = 3.2*G4Exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
      G4double Ex2 = 12*G4Exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
      fTotalXsc       = Ex1 + Ex2 + 27.5;
      fElasticXsc = 6.0 + 1.4/(( pLab - 1.4)*( pLab - 1.4) + 0.1);
    }
    else if( pLab < 2.0 ) // ns original
    {
      G4double Ex1 = 3.2*G4Exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
      G4double Ex2 = 12*G4Exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
      fTotalXsc     = Ex1 + Ex2 + 27.5;
      fElasticXsc = 3.0 + 1.36/( (logP - 0.336)*(logP - 0.336) + 0.08);    
    }
    else if( pLab < 3.5 ) // ns original
    {
      G4double Ex1 = 3.2*G4Exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
      G4double Ex2 = 12*G4Exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
      fTotalXsc       = Ex1 + Ex2 + 27.5;
      fElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
    }
    else if( pLab < 10. ) // my
    {
      fTotalXsc = 10.6 + 2.*G4Log(pE) + 25*G4Exp(-G4Log(pE)*0.43 ); 
      fElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
    }
    else //  pLab > 10 // my
    {
      fElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
    }
  }
  // pi+ n; pi- p
  else if((theParticle == thePiPlus && neutron) || 
	  (theParticle == thePiMinus && proton))  
  {
    if( pLab < 0.28 ) 
    {
      fTotalXsc = 0.288/((pLab - 0.28)*(pLab - 0.28) + 0.004);
      fElasticXsc = 1.8/((logP + 1.273)*(logP + 1.273) + 0.07);
    }
    else if( pLab < 0.395676 ) // first peak
    {
      fTotalXsc = 0.648/((pLab - 0.28)*(pLab - 0.28) + 0.009);
      fElasticXsc = 0.257/((pLab - 0.28)*(pLab - 0.28) + 0.01);
    }
    else if( pLab < 0.5 )
    {
      G4double y = G4Log(pLab/0.48);
      fTotalXsc = 26 + 110*y*y;
      fElasticXsc = 0.37*fTotalXsc;
    }
    else if( pLab < 0.65 )
    {
      G4double x = G4Log(pLab/0.48);
      fTotalXsc = 26. + 110.*x*x;
      fElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
    }
    else if( pLab < 0.72 )
    {
      fTotalXsc = 36.1 + 10*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
	24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
      fElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
    }
    else if( pLab < 0.88 )
    {
      fTotalXsc = 36.1 + 10.*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
	24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
      fElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
    }
    else if( pLab < 1.03 )
    {
      fTotalXsc = 36.1 + 10.*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
	24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
      fElasticXsc = 2.0 + 0.4/((pLab - 1.03)*(pLab - 1.03) + 0.016);
    }
    else if( pLab < 1.15 )
    {
      fTotalXsc = 36.1 + 10.*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
	24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
      fElasticXsc = 2.0 + 0.4/((pLab - 1.03)*(pLab - 1.03) + 0.016);
    }
    else if( pLab < 1.3 )
    {
      fTotalXsc = 36.1 + 10.*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
	24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
      fElasticXsc = 3. + 13./pLab;
    }
    else if( pLab < 10.) // < 3.0) // ns original
    {
      fTotalXsc = 36.1 + 0.079-4.313*logP+
	3*G4Exp(-(pLab-2.1)*(pLab-2.1)/0.4/0.4)+
	1.5*G4Exp(-(pLab-1.4)*(pLab-1.4)/0.12/0.12);
      fElasticXsc = 3. + 13./pLab; 
    }
    else   // mb 
    {
      fElasticXsc = 3. + 13./pLab;
    }
  }
  else if( (theParticle == theKMinus) && proton )   // K-p
  {
    if( pLab < pMin)
    {
      G4double psp = pLab*std::sqrt(pLab);
      fElasticXsc  = 5.2/psp;
      fTotalXsc    = 14./psp;
    }
    else if( pLab > pMax )
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc  = cofLogE*ld2 + 2.23;
      fTotalXsc    = 1.1*cofLogT*ld2 + 19.7;
    }
    else
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      G4double sp  = std::sqrt(pLab);
      G4double psp = pLab*sp;
      G4double p2  = pLab*pLab;
      G4double p4  = p2*p2;

      G4double lh  = pLab - 1.01;
      G4double hd  = lh*lh + .011;

      G4double lm  = pLab - .39;
      G4double md  = lm*lm + .000356;
      G4double lh1  = pLab - 0.78;
      G4double hd1  = lh1*lh1 + .00166;
      G4double lh2  = pLab - 1.63;
      G4double hd2  = lh2*lh2 + .007;

      // small peaks were added but commented out now
      fElasticXsc  = 5.2/psp + (1.1*cofLogE*ld2 + 2.23)/(1. - .7/sp + .075/p4) 
	+ .004/md + 0.005/hd1+ 0.01/hd2 +.15/hd; 

      fTotalXsc   = 14./psp + (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4) 
        + .006/md  + 0.01/hd1+ 0.02/hd2 + .20/hd ;
    }
  }
  else if( (theParticle == theKMinus) && neutron )  // K-n
  {
    if( pLab > pMax )
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc  = cofLogE*ld2 + 2.23;
      fTotalXsc    = 1.1*cofLogT*ld2 + 19.7;
    }
    else
    { 
      G4double lh  = pLab - 0.98;
      G4double hd  = lh*lh + .021; 
      G4double sqrLogPlab = logP*logP;

      fElasticXsc = 5.0 +  8.1*G4Exp(-logP*1.8) + 0.16*sqrLogPlab 
	- 1.3*logP + .15/hd;
      fTotalXsc = 25.2 +  0.38*sqrLogPlab - 2.9*logP + 0.30/hd;
    }
  }
  else if( (theParticle == theKPlus) && proton  )  // K+p
  {
    // VI: modified low-energy part
    if( pLab < 0.631 )
    {
      fElasticXsc = fTotalXsc = 12.03;
    }
    else if( pLab > pMax )
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc  = cofLogE*ld2 + 2.23;
      fTotalXsc    = cofLogT*ld2 + 19.2;
    }
    else
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      G4double lr  = pLab - .38;
      G4double LE  = .7/(lr*lr + .076);
      G4double sp  = std::sqrt(pLab);
      G4double p2  = pLab*pLab;
      G4double p4  = p2*p2;
      // VI: tuned elastic
      fElasticXsc  = LE + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .1/p4) 
	+ 2./((pLab - 0.8)*(pLab - 0.8) + 0.652);
      fTotalXsc    = LE + (cofLogT*ld2 + 19.5)/(1. + .46/sp + 1.6/p4) 
	+ 2.6/((pLab - 1.)*(pLab - 1.) + 0.392);
    }
  }
  else if(  (theParticle == theKPlus) && neutron) // K+n  
  {
    if( pLab < pMin )
    {
      G4double lm = pLab - 0.94;
      G4double md = lm*lm + .392;   
      fElasticXsc = 2./md;
      fTotalXsc   = 4.6/md;
    }
    else if( pLab > pMax )
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc  = cofLogE*ld2 + 2.23;
      fTotalXsc    = cofLogT*ld2 + 19.2;
    }
    else
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      G4double sp  = std::sqrt(pLab);
      G4double p2  = pLab*pLab;
      G4double p4  = p2*p2;
      G4double lm  = pLab - 0.94;
      G4double md  = lm*lm + .392; 
      fElasticXsc  = (cofLogE*ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
      fTotalXsc    = (cofLogT*ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 4.6/md;
    }
  }
  fTotalXsc   *= CLHEP::millibarn; 
  fElasticXsc *= CLHEP::millibarn;
  fElasticXsc  = std::min(fElasticXsc, fTotalXsc);

  if( proton && theParticle->GetPDGCharge() > 0. && ekin < 100*MeV)
  {
    G4double cB = G4NuclearRadii::CoulombFactor(theParticle, nucleon, ekin);
    fTotalXsc   *= cB;
    fElasticXsc *= cB; 
  }
  fInelasticXsc = std::max(fTotalXsc - fElasticXsc,0.0);
  /*
  G4cout<< "HNXsc: Ekin(GeV)= " << ekin/GeV << "; tot(mb)= " << fTotalXsc/millibarn
        <<"; el(mb)= " <<fElasticXsc/millibarn
        <<"; inel(mb)= " <<fInelasticXsc/millibarn<<"  "
        << theParticle->GetParticleName() << " + " 
        << nucleon->GetParticleName() << G4endl;
  */
  return fTotalXsc;
}

//////////////////////////////////////////////////////////////////////////////
//
// Returns kaon-nucleon cross-section based on smoothed NS 
// tuned for the Glauber-Gribov hadron model for Z>1

G4double G4HadronNucleonXsc::KaonNucleonXscGG(
                             const G4ParticleDefinition* theParticle, 
			     const G4ParticleDefinition* nucleon, G4double ekin)
{
  fTotalXsc = fElasticXsc = fInelasticXsc = 0.0;
  if(theParticle == theKMinus || theParticle == theKPlus) {
    KaonNucleonXscVG(theParticle, nucleon, ekin);

  } else if(theParticle == theK0S || theParticle == theK0L) {
    G4double stot  = KaonNucleonXscVG(theKMinus, nucleon, ekin);
    G4double sel   = fElasticXsc;
    G4double sinel = fInelasticXsc;
    stot  += KaonNucleonXscVG(theKPlus, nucleon, ekin);
    sel   += fElasticXsc;
    sinel += fInelasticXsc;
    fTotalXsc = stot*0.5; 
    fElasticXsc = sel*0.5; 
    fInelasticXsc = sinel*0.5;
  }
  return fTotalXsc;
}

//////////////////////////////////////////////////////////////////////////////
//
// Returns kaon-nucleon cross-section using NS

G4double G4HadronNucleonXsc::KaonNucleonXscNS(
                             const G4ParticleDefinition* theParticle, 
			     const G4ParticleDefinition* nucleon, G4double ekin)
{
  fTotalXsc = fElasticXsc = fInelasticXsc = 0.0;
  if(theParticle == theKMinus || theParticle == theKPlus) {
    HadronNucleonXscNS(theParticle, nucleon, ekin);

  } else if(theParticle == theK0S || theParticle == theK0L) {
    G4double stot  = HadronNucleonXscNS(theKMinus, nucleon, ekin);
    G4double sel   = fElasticXsc;
    G4double sinel = fInelasticXsc;
    stot  += HadronNucleonXscNS(theKPlus, nucleon, ekin);
    sel   += fElasticXsc;
    sinel += fInelasticXsc;
    fTotalXsc = stot*0.5; 
    fElasticXsc = sel*0.5; 
    fInelasticXsc = sinel*0.5;
  }
  return fTotalXsc;
}

//////////////////////////////////////////////////////////////////////////////
//
// Returns kaon-nucleon cross-section with smoothed NS parameterisation

G4double G4HadronNucleonXsc::KaonNucleonXscVG(
                 const G4ParticleDefinition* theParticle, 
		 const G4ParticleDefinition* nucleon, G4double ekin)
{
  G4double pM   = theParticle->GetPDGMass();
  G4double pLab = std::sqrt(ekin*(ekin + 2*pM));

  pLab *= invGeV;
  G4double logP = G4Log(pLab);

  fTotalXsc = 0.0;

  G4bool proton = (nucleon == theProton);
  G4bool neutron = (nucleon == theNeutron);

  if( (theParticle == theKMinus) && proton )   // K-p
  {
    if( pLab < pMin)
    {
      G4double psp = pLab*std::sqrt(pLab);
      fElasticXsc  = 5.2/psp;
      fTotalXsc    = 14./psp;
    }
    else if( pLab > pMax )
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc  = cofLogE*ld2 + 2.23;
      fTotalXsc    = 1.1*cofLogT*ld2 + 19.7;
    }
    else
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      G4double sp  = std::sqrt(pLab);
      G4double psp = pLab*sp;
      G4double p2  = pLab*pLab;
      G4double p4  = p2*p2;

      G4double lh  = pLab - 1.01;
      G4double hd  = lh*lh + .011;
      fElasticXsc  = 5.2/psp + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .075/p4) + .15/hd;
      fTotalXsc    = 14./psp + (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4) + .60/hd;
    }
  }
  else if( (theParticle == theKMinus) && neutron )  // K-n
  {
    if( pLab > pMax )
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc  = cofLogE*ld2 + 2.23;
      fTotalXsc    = 1.1*cofLogT*ld2 + 19.7;
    }
    else
    { 
      G4double lh  = pLab - 0.98;
      G4double hd  = lh*lh + .045;    // vg version
      G4double sqrLogPlab = logP*logP;

      fElasticXsc = 5.0 +  8.1*G4Exp(-logP*1.8) + 0.16*sqrLogPlab 
	- 1.3*logP + .15/hd;
      fTotalXsc = 25.2 +  0.38*sqrLogPlab - 2.9*logP + 0.60/hd;  // vg version
    }
  }
  else if( (theParticle == theKPlus) && proton  )  // K+p
  {
    // VI: modified low-energy part
    if( pLab < 0.631 )
    {
      fElasticXsc = fTotalXsc = 12.03;
    }
    else if( pLab > pMax )
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc  = cofLogE*ld2 + 2.23;
      fTotalXsc    = cofLogT*ld2 + 19.2;
    }
    else
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      G4double lr  = pLab - .38;
      G4double LE  = .7/(lr*lr + .076);
      G4double sp  = std::sqrt(pLab);
      G4double p2  = pLab*pLab;
      G4double p4  = p2*p2;
      // VI: tuned elastic
      fElasticXsc  = LE + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .1/p4) 
	+ 2./((pLab - 0.8)*(pLab - 0.8) + 0.652);
      fTotalXsc    = LE + (cofLogT*ld2 + 19.5)/(1. + .46/sp + 1.6/p4) 
	+ 2.6/((pLab - 1.)*(pLab - 1.) + 0.392);
    }
  }
  else if(  (theParticle == theKPlus) && neutron) // K+n  
  {
    if( pLab < pMin )
    {
      G4double lm = pLab - 0.94;
      G4double md = lm*lm + .392;   
      fElasticXsc = 2./md;
      fTotalXsc   = 4.6/md;
    }
    else if( pLab > pMax )
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc  = cofLogE*ld2 + 2.23;
      fTotalXsc    = cofLogT*ld2 + 19.2;
    }
    else
    {
      G4double ld  = logP - minLogP;
      G4double ld2 = ld*ld;
      G4double sp  = std::sqrt(pLab);
      G4double p2  = pLab*pLab;
      G4double p4  = p2*p2;
      G4double lm  = pLab - 0.94;
      G4double md  = lm*lm + .392; 
      fElasticXsc  = (cofLogE*ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
      fTotalXsc    = (cofLogT*ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 4.6/md;
    }
  }

  fTotalXsc   *= CLHEP::millibarn; 
  fElasticXsc *= CLHEP::millibarn;

  if( proton && theParticle->GetPDGCharge() > 0. )
  {
    G4double cB = G4NuclearRadii::CoulombFactor(theParticle, nucleon, ekin);
    fTotalXsc   *= cB;
    fElasticXsc *= cB; 
  }
  fElasticXsc = std::min(fElasticXsc, fTotalXsc);
  fInelasticXsc = std::max(fTotalXsc - fElasticXsc,0.0);
  /*
  G4cout << "HNXscVG: E= " << ekin << " " << theParticle->GetParticleName() 
	 << " P: " << proton << " xtot(b)= " << fTotalXsc/barn
	 << " xel(b)= " <<  fElasticXsc/barn << " xinel(b)= " << fInelasticXsc/barn
	 << G4endl;
  */
  return fTotalXsc;
}

//////////////////////////////////////////////////////////////////////////////
//
// Returns hyperon-nucleon cross-section using NS x-section for protons

G4double G4HadronNucleonXsc::HyperonNucleonXscNS(
         const G4ParticleDefinition* theParticle, 
	 const G4ParticleDefinition* nucleon, G4double ekin)
{
  G4double coeff = 1.0;
  
  static const G4double  lBarCof1S  = 0.88;
  static const G4double  lBarCof2S  = 0.76;
  static const G4double  lBarCof3S  = 0.64;
  static const G4double  lBarCof1C  = 0.784378;
  static const G4double  lBarCofSC  = 0.664378;
  static const G4double  lBarCof2SC = 0.544378;
  static const G4double  lBarCof1B  = 0.740659;
  static const G4double  lBarCofSB  = 0.620659;
  static const G4double  lBarCof2SB = 0.500659;
  
  if( theParticle == theL || theParticle == theSPlus || 
      theParticle == theSMinus || theParticle == theS0 ||
      theParticle == theAntiL || theParticle == theASPlus || 
      theParticle == theASMinus || theParticle == theAS0  )
  {
    coeff = lBarCof1S;

  } else if( theParticle == theXiMinus || theParticle == theXi0 || 
             theParticle == theAXiMinus || theParticle == theAXi0  )
  {
    coeff = lBarCof2S;
  }
  else if( theParticle == theOmega || theParticle == theAOmega)
  {
    coeff = lBarCof3S;
  }
  else if( theParticle == theLambdaCPlus || theParticle == theALambdaCPlus ||
	   theParticle == theSigmaCPlus || theParticle == theASigmaCPlus ||
	   theParticle == theSigmacPP || theParticle == theASigmacPP     ||
	   theParticle == theSigmaC0 || theParticle == theASigmaC0
	 )
  {
    coeff = lBarCof1C;
  }
  else if( theParticle == theOmegaC0 || theParticle == theAOmegaC0 )
  {
    coeff = lBarCof2SC;
  }
  else if( theParticle == theXiCPlus || theParticle == theXiC0 || 
             theParticle == theAXiCPlus || theParticle == theAXiC0)
  {
    coeff = lBarCofSC;
  }
  else if( theParticle == theLambdaB || theParticle == theALambdaB ||
	   theParticle == theSigmaBPlus || theParticle == theASigmaBPlus ||
	   theParticle == theSigmaBMinus || theParticle == theASigmaBMinus     ||
	   theParticle == theSigmaB0 || theParticle == theASigmaB0
	 )
  {
    coeff = lBarCof1B;
  }
  else if( theParticle == theOmegaBMinus || theParticle == theAOmegaBMinus)
  {
    coeff = lBarCof2SB;
  }
  else if( theParticle == theXiBMinus || theParticle == theXiB0 || 
             theParticle == theAXiBMinus || theParticle == theAXiB0)
  {
    coeff = lBarCofSB;
  } 
  fTotalXsc = coeff*HadronNucleonXscNS( theProton, nucleon, ekin);
  fInelasticXsc *= coeff;
  fElasticXsc *= coeff;
  
  return fTotalXsc;
}

//////////////////////////////////////////////////////////////////////////////
//
// Returns hyperon-nucleon cross-section using NS x-section for protons

G4double G4HadronNucleonXsc::SCBMesonNucleonXscNS( 
         const G4ParticleDefinition* theParticle, 
         const G4ParticleDefinition* nucleon, G4double ekin )
{
  G4double coeff(1.0);
  // static const G4double lMesCof1S = 0.82; // Kp/piP
  static const G4double llMesCof1C = 0.676568;
  static const G4double llMesCof1B = 0.610989;
  static const G4double llMesCof2C = 0.353135;
  static const G4double llMesCof2B = 0.221978;
  static const G4double llMesCofSC = 0.496568;
  static const G4double llMesCofSB = 0.430989;
  static const G4double llMesCofCB = 0.287557;
  static const G4double llMesCofEtaP = 0.88;
  static const G4double llMesCofEta = 0.76;

  if( theParticle == theBMeson0 || theParticle == theABMeson0 ||
      theParticle == theBMesonPlus || theParticle == theBMesonMinus )
  {
    coeff = llMesCof1B;
  }
  else if(theParticle == theDMeson0 || theParticle == theADMeson0 ||
      theParticle == theDMesonPlus || theParticle == theDMesonMinus )
  {
    coeff = llMesCof1C;
  }
  else if(theParticle == theBsMeson0 || theParticle == theABsMeson0 )
  {
    coeff = llMesCofSB;
  }
  else if(theParticle == theBcMesonPlus || theParticle == theBcMesonMinus )
  {
    coeff = llMesCofCB;
  }
  else if(theParticle == theDsMesonPlus || theParticle == theDsMesonMinus )
  {
    coeff = llMesCofSC;
  }
  else if(theParticle == theBMesonPlus || theParticle == theBMesonMinus )
  {
    coeff = llMesCof1B;
  }
  else if(theParticle == theDMesonPlus || theParticle == theDMesonMinus )
  {
    coeff = llMesCof1C;
  }
  else if(theParticle == theEtaC || theParticle == theJPsi )
  {
    coeff = llMesCof2C;
  }
  else if(theParticle == theUpsilon )
  {
    coeff = llMesCof2B;
  }
  else if(theParticle == theEta )
  {
    coeff = llMesCofEta;
  }
  else if(theParticle == theEtaPrime )
  {
    coeff = llMesCofEtaP;
  }
  fTotalXsc = coeff*HadronNucleonXscNS( thePiPlus, nucleon, ekin);  
  fElasticXsc *= coeff;
  fInelasticXsc *= coeff;
  return fTotalXsc;
}
////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon cross-section based on V. Uzjinsky parametrisation of
// data from G4FTFCrossSection class

G4double G4HadronNucleonXsc::HadronNucleonXscVU(
                             const G4ParticleDefinition* theParticle, 
			     const G4ParticleDefinition* nucleon, G4double ekin)
{
  G4int PDGcode = theParticle->GetPDGEncoding();
  G4int absPDGcode = std::abs(PDGcode);
  G4double mass = theParticle->GetPDGMass();
  G4double Elab = ekin + mass;              
  // (s - 2*0.88*GeV*GeV)/(2*0.939*GeV)/GeV;
  G4double Plab = std::sqrt(ekin*(ekin + 2.*mass));            

  Elab *= invGeV;
  Plab *= invGeV;

  G4double logPlab = G4Log( Plab );
  G4double sqrLogPlab = logPlab * logPlab;

  G4bool proton = (nucleon == theProton);
  G4bool neutron = (nucleon == theNeutron);
   
  if( absPDGcode > 1000)  //------Projectile is baryon -
  {
    if(proton)
    {
      fTotalXsc   = 48.0 + 0.522*sqrLogPlab - 4.51*logPlab;
      fElasticXsc = 11.9 + 26.9*G4Exp(-logPlab*1.21) + 0.169*sqrLogPlab - 1.85*logPlab;
    }
    if(neutron)
    {    
      fTotalXsc   = 47.3  + 0.513*sqrLogPlab - 4.27*logPlab;
      fElasticXsc = 11.9 + 26.9*G4Exp(-logPlab*1.21) + 0.169*sqrLogPlab - 1.85*logPlab;
    }
  }
  else if( PDGcode ==  211)  //------Projectile is PionPlus ----
  {
    if(proton)
    {
      fTotalXsc  = 16.4 + 19.3 *G4Exp(-logPlab*0.42) + 0.19 *sqrLogPlab - 0.0 *logPlab;
      fElasticXsc =  0.0 + 11.4*G4Exp(-logPlab*0.40) + 0.079*sqrLogPlab - 0.0 *logPlab;
    }
    if(neutron)
    {    
      fTotalXsc   =  33.0 + 14.0 *G4Exp(-logPlab*1.36) + 0.456*sqrLogPlab - 4.03*logPlab;
      fElasticXsc = 1.76 + 11.2*G4Exp(-logPlab*0.64) + 0.043*sqrLogPlab - 0.0 *logPlab;
    }
  }
  else if( PDGcode == -211)  //------Projectile is PionMinus ----
  {
    if(proton)
    {
      fTotalXsc   = 33.0 + 14.0 *G4Exp(-logPlab*1.36) + 0.456*sqrLogPlab - 4.03*logPlab;
      fElasticXsc = 1.76 + 11.2*G4Exp(-logPlab*0.64) + 0.043*sqrLogPlab - 0.0 *logPlab;
    }
    if(neutron)
    {    
      fTotalXsc   = 16.4 + 19.3 *G4Exp(-logPlab*0.42) + 0.19 *sqrLogPlab - 0.0 *logPlab;
      fElasticXsc =  0.0 + 11.4*G4Exp(-logPlab*0.40) + 0.079*sqrLogPlab - 0.0 *logPlab;
    }
  }
  else if( PDGcode ==  111)  //------Projectile is PionZero  --
  {
    if(proton)
    {
      fTotalXsc = (16.4 + 19.3 *G4Exp(-logPlab*0.42) + 0.19 *sqrLogPlab - 0.0 *logPlab +   //Pi+
		   33.0 + 14.0 *G4Exp(-logPlab*1.36) + 0.456*sqrLogPlab - 4.03*logPlab)/2; //Pi-

      fElasticXsc = (0.0 + 11.4*G4Exp(-logPlab*0.40) + 0.079*sqrLogPlab - 0.0 *logPlab +   //Pi+
		     1.76 + 11.2*G4Exp(-logPlab*0.64) + 0.043*sqrLogPlab - 0.0 *logPlab)/2;//Pi-

    }
    if(neutron)
    {    
      fTotalXsc = (33.0 + 14.0 *G4Exp(-logPlab*1.36) + 0.456*sqrLogPlab - 4.03*logPlab +   //Pi+
		   16.4 + 19.3 *G4Exp(-logPlab*0.42) + 0.19 *sqrLogPlab - 0.0 *logPlab)/2; //Pi-
      fElasticXsc = (1.76 + 11.2*G4Exp(-logPlab*0.64) + 0.043*sqrLogPlab - 0.0 *logPlab +  //Pi+
		     0.0  + 11.4*G4Exp(-logPlab*0.40) + 0.079*sqrLogPlab - 0.0 *logPlab)/2;//Pi-
    }
  }
  else if( PDGcode == 321 )    //------Projectile is KaonPlus --
  {
    if(proton)
    {
      fTotalXsc   = 18.1 +  0.26 *sqrLogPlab - 1.0 *logPlab;
      fElasticXsc =  5.0 +  8.1*G4Exp(-logPlab*1.8 ) + 0.16 *sqrLogPlab - 1.3 *logPlab;
    }
    if(neutron)
    {    
      fTotalXsc   = 18.7  + 0.21 *sqrLogPlab - 0.89*logPlab;
      fElasticXsc =  7.3  + 0.29 *sqrLogPlab - 2.4 *logPlab;
    }
  }
  else if( PDGcode ==-321 )  //------Projectile is KaonMinus ----
  {
    if(proton)
    {
      fTotalXsc   = 32.1 + 0.66*sqrLogPlab - 5.6*logPlab;
      fElasticXsc =  7.3 + 0.29*sqrLogPlab - 2.4*logPlab;
    }
    if(neutron)
    {    
      fTotalXsc   = 25.2 + 0.38*sqrLogPlab - 2.9*logPlab;
      fElasticXsc =  5.0 +  8.1*G4Exp(-logPlab*1.8 ) + 0.16*sqrLogPlab - 1.3*logPlab;
    }
  }
  else if( PDGcode == 311 )  //------Projectile is KaonZero -----
  {
    if(proton)
    {
      fTotalXsc   = ( 18.1 + 0.26 *sqrLogPlab - 1.0 *logPlab +   //K+
                      32.1 + 0.66 *sqrLogPlab - 5.6 *logPlab)/2; //K-
      fElasticXsc = (  5.0 +  8.1*G4Exp(-logPlab*1.8 ) + 0.16 *sqrLogPlab - 1.3 *logPlab + //K+
                         7.3 + 0.29 *sqrLogPlab - 2.4 *logPlab)/2; //K-
    }
    if(neutron)
    {    
      fTotalXsc   = ( 18.7 + 0.21 *sqrLogPlab - 0.89*logPlab +   //K+
                      25.2 + 0.38 *sqrLogPlab - 2.9 *logPlab)/2; //K-
      fElasticXsc = (  7.3 + 0.29 *sqrLogPlab - 2.4 *logPlab +   //K+
                       5.0 + 8.1*G4Exp(-logPlab*1.8 ) + 0.16 *sqrLogPlab - 1.3 *logPlab)/2; //K-
    }
  }
  else  //------Projectile is undefined, Nucleon assumed
  {
    if(proton)
    {
      fTotalXsc   = 48.0 + 0.522*sqrLogPlab - 4.51*logPlab;
      fElasticXsc = 11.9 + 26.9*G4Exp(-logPlab*1.21) + 0.169*sqrLogPlab - 1.85*logPlab;
    }
    if(neutron)
    {    
      fTotalXsc   = 47.3 + 0.513*sqrLogPlab - 4.27*logPlab;
      fElasticXsc = 11.9 + 26.9*G4Exp(-logPlab*1.21) + 0.169*sqrLogPlab - 1.85*logPlab;
    }
  }

  fTotalXsc   *= CLHEP::millibarn;
  fElasticXsc *= CLHEP::millibarn;
  fElasticXsc = std::min(fElasticXsc, fTotalXsc);
  fInelasticXsc = fTotalXsc - fElasticXsc;

  return fTotalXsc;    
}

//////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to differnt parametrisations:
// [2] E. Levin, hep-ph/9710546
// [3] U. Dersch, et al, hep-ex/9910052
// [4] M.J. Longo, et al, Phys.Rev.Lett. 33 (1974) 725 

G4double G4HadronNucleonXsc::HadronNucleonXscEL(
                             const G4ParticleDefinition* theParticle, 
			     const G4ParticleDefinition*, G4double ekin)
{
  G4double xsection(0.);
  static const G4double targ_mass = 
    0.5*(CLHEP::proton_mass_c2 + CLHEP::neutron_mass_c2);
  G4double sMand = 
    CalcMandelstamS(ekin, theParticle->GetPDGMass(), targ_mass)*invGeV2;

  G4double x1 = G4Exp(G4Log(sMand)*0.0808);
  G4double x2 = G4Exp(G4Log(-sMand)*0.4525);

  if(theParticle == theGamma) 
  {
    xsection = 0.0677*x1 + 0.129*x2;
  } 
  else if(theParticle == theNeutron)  
  {
    xsection = 21.70*x1 + 56.08*x2;
  } 
  else if(theParticle == theProton) 
  {
    xsection = 21.70*x1 + 56.08*x2;
  } 
  else if(theParticle == theAProton) 
  {
    xsection = 21.70*x1 + 98.39*x2;
  } 
  else if(theParticle == thePiPlus) 
  {
    xsection = 13.63*x1 + 27.56*x2;
  } 
  else if(theParticle == thePiMinus) 
  {
    xsection = 13.63*x1 + 36.02*x2;
  } 
  else if(theParticle == theKPlus) 
  {
    xsection = 11.82*x1 + 8.15*x2;
  } 
  else if(theParticle == theKMinus) 
  {
    xsection = 11.82*x1 + 26.36*x2;
  }
  else if(theParticle == theK0S || theParticle == theK0L) 
  {
    xsection = 11.82*x1 + 17.25*x2;
  }
  else  
  {
    xsection = 21.70*x1 + 56.08*x2;
  }
  fTotalXsc     = xsection*CLHEP::millibarn;
  fInelasticXsc = 0.83*fTotalXsc;
  fElasticXsc   = fTotalXsc - fInelasticXsc;
  return fTotalXsc;
}

///////////////////////////////////////////////////////////////////////////////

G4double 
G4HadronNucleonXsc::CoulombBarrier(const G4ParticleDefinition* theParticle, 
				   const G4ParticleDefinition* nucleon, 
				   G4double ekin)
{
  G4double tR = 0.895*CLHEP::fermi;
  G4double pR = 0.5*CLHEP::fermi;

  if     ( theParticle == theProton ) pR = 0.895*fermi;
  else if( theParticle == thePiPlus ) pR = 0.663*fermi;
  else if( theParticle == theKPlus )  pR = 0.340*fermi;

  G4double pZ = theParticle->GetPDGCharge();
  G4double tZ = nucleon->GetPDGCharge();

  G4double pM = theParticle->GetPDGMass(); 
  G4double tM = nucleon->GetPDGMass();

  G4double pElab = ekin + pM;

  G4double totEcm  = std::sqrt(pM*pM + tM*tM + 2.*pElab*tM);

  G4double totTcm  = totEcm - pM -tM;

  G4double bC = fine_structure_const*hbarc*pZ*tZ/(2.*(pR + tR));

  //G4cout<<"##CB: ekin = "<<ekin<<"; pElab= " << pElab
  //	<<"; bC = "<<bC<<"; Tcm = "<< totTcm <<G4endl;

  G4double ratio = (totTcm > bC) ? 1. - bC/totTcm : 0.0;

  // G4cout <<"ratio = "<<ratio<<G4endl;
  return ratio;
}

//////////////////////////////////////////////////////////////////////////////

