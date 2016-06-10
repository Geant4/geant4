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

#include "G4HadronNucleonXsc.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadTmpUtil.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"

G4HadronNucleonXsc::G4HadronNucleonXsc() 
: 
// fUpperLimit( 10000 * GeV ),
  fLowerLimit( 0.03 * MeV ),
  fTotalXsc(0.0), fElasticXsc(0.0), fInelasticXsc(0.0)
// , fHadronNucleonXsc(0.0)
{
  theGamma    = G4Gamma::Gamma();
  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theAProton  = G4AntiProton::AntiProton();
  theANeutron = G4AntiNeutron::AntiNeutron();
  thePiPlus   = G4PionPlus::PionPlus();
  thePiMinus  = G4PionMinus::PionMinus();
  thePiZero   = G4PionZero::PionZero();
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
  theD        = G4Deuteron::Deuteron();
  theT        = G4Triton::Triton();
  theA        = G4Alpha::Alpha();
  theHe3      = G4He3::He3();

  // InitialiseKaonNucleonTotXsc();
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

G4bool 
G4HadronNucleonXsc::IsApplicable(const G4DynamicParticle* aDP, 
                                 const G4Element* anElement)
{
  G4int Z = G4lrint(anElement->GetZ());
  G4int A = G4lrint(anElement->GetN());
  return IsIsoApplicable(aDP, Z, A);
} 

////////////////////////////////////////////////////////////////////////////////////////
//

G4bool 
G4HadronNucleonXsc::IsIsoApplicable(const G4DynamicParticle* aDP, 
                                    G4int Z, G4int)
{
  G4bool applicable = false;
  // G4int baryonNumber     = aDP->GetDefinition()->GetBaryonNumber();
  G4double kineticEnergy = aDP->GetKineticEnergy();

  const G4ParticleDefinition* theParticle = aDP->GetDefinition();
 
  if ( ( kineticEnergy  >= fLowerLimit &&
         Z > 1 &&      // >=  He
       ( theParticle == theAProton   ||
         theParticle == theGamma     ||
         theParticle == theKPlus     ||
         theParticle == theKMinus    || 
         theParticle == theSMinus)      )    ||  

       ( kineticEnergy  >= 0.1*fLowerLimit &&
         Z > 1 &&      // >=  He
       ( theParticle == theProton    ||
         theParticle == theNeutron   ||   
         theParticle == thePiPlus    ||
         theParticle == thePiMinus       ) )    ) applicable = true;

  return applicable;
}


/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to differnt parametrisations:
// [2] E. Levin, hep-ph/9710546
// [3] U. Dersch, et al, hep-ex/9910052
// [4] M.J. Longo, et al, Phys.Rev.Lett. 33 (1974) 725 

G4double 
G4HadronNucleonXsc::GetHadronNucleonXscEL(const G4DynamicParticle* aParticle, 
                                          const G4ParticleDefinition* nucleon )
{
  G4double xsection;


  G4double targ_mass = 0.939*GeV;  // ~mean neutron and proton ???

  G4double proj_mass     = aParticle->GetMass();
  G4double proj_momentum = aParticle->GetMomentum().mag();
  G4double sMand = CalcMandelstamS ( proj_mass , targ_mass , proj_momentum );

  sMand /= GeV*GeV;  // in GeV for parametrisation
  proj_momentum /= GeV;

  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();

  G4bool pORn = (nucleon == theProton || nucleon == theNeutron  );  
  

  if(theParticle == theGamma && pORn ) 
  {
    xsection = (0.0677*std::pow(sMand,0.0808) + 0.129*std::pow(sMand,-0.4525));
  } 
  else if(theParticle == theNeutron && pORn ) // as proton ??? 
  {
    xsection = (21.70*std::pow(sMand,0.0808) + 56.08*std::pow(sMand,-0.4525));
  } 
  else if(theParticle == theProton && pORn ) 
  {
    xsection = (21.70*std::pow(sMand,0.0808) + 56.08*std::pow(sMand,-0.4525));

    // xsection = At*( 49.51*std::pow(sMand,-0.097) + 0.314*G4Log(sMand)*G4Log(sMand) );
    // xsection = At*( 38.4 + 0.85*std::abs(std::pow(log(sMand),1.47)) );
  } 
  else if(theParticle == theAProton && pORn ) 
  {
    xsection = ( 21.70*std::pow(sMand,0.0808) + 98.39*std::pow(sMand,-0.4525));
  } 
  else if(theParticle == thePiPlus && pORn ) 
  {
    xsection = (13.63*std::pow(sMand,0.0808) + 27.56*std::pow(sMand,-0.4525));
  } 
  else if(theParticle == thePiMinus && pORn ) 
  {
    // xsection = At*( 55.2*std::pow(sMand,-0.255) + 0.346*G4Log(sMand)*G4Log(sMand) );
    xsection = (13.63*std::pow(sMand,0.0808) + 36.02*std::pow(sMand,-0.4525));
  } 
  else if(theParticle == theKPlus && pORn ) 
  {
    xsection = (11.82*std::pow(sMand,0.0808) + 8.15*std::pow(sMand,-0.4525));
  } 
  else if(theParticle == theKMinus && pORn ) 
  {
    xsection = (11.82*std::pow(sMand,0.0808) + 26.36*std::pow(sMand,-0.4525));
  }
  else  // as proton ??? 
  {
    xsection = (21.70*std::pow(sMand,0.0808) + 56.08*std::pow(sMand,-0.4525));
  }
  xsection *= millibarn;

  fTotalXsc     = xsection;
  fInelasticXsc = 0.83*xsection;
  fElasticXsc   = fTotalXsc - fInelasticXsc;
  if (fElasticXsc < 0.)fElasticXsc = 0.;
 
  return xsection;
}


/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to PDG parametrisation (2005):
// http://pdg.lbl.gov/2006/reviews/hadronicrpp.pdf
//  At = number of nucleons,  Zt = number of protons 

G4double 
G4HadronNucleonXsc::GetHadronNucleonXscPDG(const G4DynamicParticle* aParticle, 
                                           const G4ParticleDefinition* nucleon )
{
  G4double xsection(0);
  G4int Zt=1, Nt=1, At=1;

   G4double targ_mass = 0.939*GeV;  // ~mean neutron and proton ???

  G4double proj_mass     = aParticle->GetMass(); 
  G4double proj_momentum = aParticle->GetMomentum().mag();

  G4double sMand = CalcMandelstamS ( proj_mass , targ_mass , proj_momentum );

  sMand         /= GeV*GeV;  // in GeV for parametrisation

  // General PDG fit constants

  G4double s0   = 5.38*5.38; // in Gev^2
  G4double eta1 = 0.458;
  G4double eta2 = 0.458;
  G4double B    = 0.308;

  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();

  G4bool pORn = (nucleon == theProton || nucleon == theNeutron  );  
  G4bool proton = (nucleon == theProton);
  G4bool neutron = (nucleon == theNeutron);
  
  if(theParticle == theNeutron) // proton-neutron fit 
  {
    if ( proton )
    {
      xsection = Zt*( 35.80 + B*std::pow(G4Log(sMand/s0),2.) 
		 + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));// on p
    }
    if ( neutron )
    {
      xsection  = Nt*( 35.45 + B*std::pow(G4Log(sMand/s0),2.) 
		      + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2)); // on n pp for nn
    }
  } 
  else if(theParticle == theProton) 
  {
    if ( proton )
    {      
      xsection  = Zt*( 35.45 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = Nt*( 35.80 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
    }
  } 
  else if(theParticle == theAProton) 
  {
    if ( proton )
    {      
      xsection  = Zt*( 35.45 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) + 33.34*std::pow(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = Nt*( 35.80 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) + 30.*std::pow(sMand,-eta2));
    }
  } 
  else if(theParticle == thePiPlus && pORn ) 
  {
    xsection  = At*( 20.86 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 19.24*std::pow(sMand,-eta1) - 6.03*std::pow(sMand,-eta2));
  } 
  else if(theParticle == thePiMinus && pORn ) 
  {
    xsection  = At*( 20.86 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 19.24*std::pow(sMand,-eta1) + 6.03*std::pow(sMand,-eta2));
  } 
  else if(theParticle == theKPlus) 
  {
    if ( proton )
    {      
      xsection  = Zt*( 17.91 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 7.14*std::pow(sMand,-eta1) - 13.45*std::pow(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = Nt*( 17.87 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 5.17*std::pow(sMand,-eta1) - 7.23*std::pow(sMand,-eta2));
    }
  } 
  else if(theParticle == theKMinus) 
  {
    if ( proton )
    {      
      xsection  = Zt*( 17.91 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 7.14*std::pow(sMand,-eta1) + 13.45*std::pow(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = Nt*( 17.87 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 5.17*std::pow(sMand,-eta1) + 7.23*std::pow(sMand,-eta2) );
    }
  }
  else if(theParticle == theSMinus && pORn ) 
  {
    xsection  = At*( 35.20 + B*std::pow(G4Log(sMand/s0),2.) 
                          - 199.*std::pow(sMand,-eta1) + 264.*std::pow(sMand,-eta2) );
  } 
  else if(theParticle == theGamma && pORn ) // modify later on
  {
    xsection  = At*( 0.0 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 0.032*std::pow(sMand,-eta1) - 0.0*std::pow(sMand,-eta2) );
   
  } 
  else  // as proton ??? 
  {
    if ( proton )
    {      
      xsection  = Zt*( 35.45 + B*std::pow(G4Log(sMand/s0),2.) 
                       + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2) );
    }
    if ( neutron )
    {
      xsection = Nt*( 35.80 + B*std::pow(G4Log(sMand/s0),2.) 
                      + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
    }
  } 
  xsection *= millibarn; // parametrised in mb

  fTotalXsc     = xsection;
  fInelasticXsc = 0.75*xsection;
  fElasticXsc   = fTotalXsc - fInelasticXsc;
  if (fElasticXsc < 0.) fElasticXsc = 0.;

  return xsection;
}


/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon cross-section based on N. Starkov parametrisation of
// data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database

G4double 
G4HadronNucleonXsc::GetHadronNucleonXscNS(const G4DynamicParticle* aParticle, 
                                          const G4ParticleDefinition* nucleon  )
{
  G4double xsection(0); 
  
  G4double A0, B0;
  G4double hpXsc(0);
  G4double hnXsc(0);


  G4double tM = 0.939*GeV;  // ~mean neutron and proton ???

  G4double pM   = aParticle->GetMass();
  G4double pE   = aParticle->GetTotalEnergy(); 
  G4double pLab = aParticle->GetMomentum().mag();

  G4double sMand = CalcMandelstamS ( pM , tM , pLab );

  sMand         /= GeV*GeV;  // in GeV for parametrisation
  pLab /= GeV;
  pE   /= GeV;
  pM     /= GeV;

  G4double logP = G4Log(pLab);


  // General PDG fit constants

  G4double s0   = 5.38*5.38; // in Gev^2
  G4double eta1 = 0.458;
  G4double eta2 = 0.458;
  G4double B    = 0.308;
  G4double minLogP = 3.5;       // min of (lnP-minLogP)^2 
  G4double cofLogE = .0557;     // elastic (lnP-minLogP)^2 
  G4double cofLogT = .3;        // total (lnP-minLogP)^2 
  G4double pMin = .1;        // fast LE calculation 
  G4double pMax = 1000.;     // fast HE calculation 


  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();

  G4bool pORn = (nucleon == theProton || nucleon == theNeutron  );  
  G4bool proton = (nucleon == theProton);
  G4bool neutron = (nucleon == theNeutron);

  if( theParticle == theNeutron && pORn ) 
  {
    if( pLab >= 373.)
    {
      xsection =  GetHadronNucleonXscPDG(aParticle, nucleon)/millibarn;

      fElasticXsc = 6.5 + 0.308*std::pow(G4Log(sMand/400.),1.65) + 9.19*std::pow(sMand,-0.458);
    
      fTotalXsc = xsection;
    
    }
    else if( pLab >= 100.)
    {
      B0 = 7.5;
      A0 = 100. - B0*G4Log(3.0e7);

      xsection = A0 + B0*G4Log(pE) - 11
	  // + 103*std::pow(2*0.93827*pE + pM*pM+0.93827*0.93827,-0.165);        //  mb
                  + 103*std::pow(sMand,-0.165);        //  mb

      fElasticXsc = 5.53 + 0.308*std::pow(G4Log(sMand/28.9),1.1) + 9.19*std::pow(sMand,-0.458);
      
      fTotalXsc = xsection;
    }
    else if( pLab >= 10.)
    {
        B0 = 7.5;
        A0 = 100. - B0*G4Log(3.0e7);

        xsection = A0 + B0*G4Log(pE) - 11
                  + 103*std::pow(2*0.93827*pE + pM*pM+
                     0.93827*0.93827,-0.165);        //  mb      
      fTotalXsc = xsection;
      fElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
    }
    else  // pLab < 10 GeV/c
    {
      if( neutron )      // nn to be pp
      {
        if( pLab < 0.4 )
        {
          hnXsc = 23 + 50*( std::pow( G4Log(0.73/pLab), 3.5 ) );
          fElasticXsc = hnXsc;
        }
        else if( pLab < 0.73 )
        {
          hnXsc = 23 + 50*( std::pow( G4Log(0.73/pLab), 3.5 ) );
          fElasticXsc = hnXsc; 
        }
        else if( pLab < 1.05  )
        {
          hnXsc = 23 + 40*(G4Log(pLab/0.73))*
                         (G4Log(pLab/0.73));
          fElasticXsc = 23 + 20*(G4Log(pLab/0.73))*
                         (G4Log(pLab/0.73));
        }
        else    // 1.05 - 10 GeV/c
        {
          hnXsc = 39.0+75*(pLab - 1.2)/(std::pow(pLab,3.0) + 0.15);

          fElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
        }
        fTotalXsc = hnXsc;
      }
      if( proton )   // pn to be np
      {
        if( pLab < 0.02 )
        {
          hpXsc = 4100+30*std::pow(G4Log(1.3/pLab),3.6); // was as pLab < 0.8
	  fElasticXsc = hpXsc;
        }      
        else if( pLab < 0.8 )
        {
          hpXsc = 33+30*std::pow(G4Log(pLab/1.3),4.0);
	  fElasticXsc = hpXsc;
        }      
        else if( pLab < 1.05 )
        {
          hpXsc = 33+30*std::pow(G4Log(pLab/0.95),2.0);
          fElasticXsc =  6 + 52/( G4Log(0.511/pLab)*G4Log(0.511/pLab) + 1.6 );
        }
        else if( pLab < 1.4 )
        {
          hpXsc = 33+30*std::pow(G4Log(pLab/0.95),2.0);
          fElasticXsc =  6 + 52/( G4Log(0.511/pLab)*G4Log(0.511/pLab) + 1.6 );
        }
        else    // 1.4 < pLab < 10.  )
        {
          hpXsc = 33.3 + 20.8*(std::pow(pLab,2.0) - 1.35)/(std::pow(pLab,2.50) + 0.95);
          
          fElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
        }
        fTotalXsc = hpXsc;
      }
    }
  } 
  else if( theParticle == theProton && pORn ) ////// proton //////////////////////////////////////////////
  {
    if( pLab >= 373.) // pdg due to TOTEM data
    {
      xsection =  GetHadronNucleonXscPDG(aParticle, nucleon)/millibarn;

      fElasticXsc = 6.5 + 0.308*std::pow(G4Log(sMand/400.),1.65) + 9.19*std::pow(sMand,-0.458);
     
      fTotalXsc = xsection;
    }
    else if( pLab >= 100.)
    {
      B0 = 7.5;
      A0 = 100. - B0*G4Log(3.0e7);

      xsection = A0 + B0*G4Log(pE) - 11 + 103*std::pow(sMand,-0.165);        //  mb

      fElasticXsc = 5.53 + 0.308*std::pow(G4Log(sMand/28.9),1.1) + 9.19*std::pow(sMand,-0.458);
      
      fTotalXsc = xsection;
    }
    else if( pLab >= 10.)
    {
      B0 = 7.5;
      A0 = 100. - B0*G4Log(3.0e7);

      xsection = A0 + B0*G4Log(pE) - 11 + 103*std::pow(sMand,-0.165);        //  mb

      fElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
      
      fTotalXsc = xsection;
    }
    else
    {
      // pp

      if( proton )
      {
        if( pLab < 0.4 )
        {
          hpXsc = 23 + 50*( std::pow( G4Log(0.73/pLab), 3.5 ) );
          fElasticXsc = hpXsc;
        }
        else if( pLab < 0.73 )
        {
          hpXsc = 23 + 50*( std::pow( G4Log(0.73/pLab), 3.5 ) );
          fElasticXsc = hpXsc; 
        }
        else if( pLab < 1.05  )
        {
          hpXsc = 23 + 40*(G4Log(pLab/0.73))*
                         (G4Log(pLab/0.73));
          fElasticXsc = 23 + 20*(G4Log(pLab/0.73))*
                         (G4Log(pLab/0.73));
        }
        else    // 1.05 - 10 GeV/c
        {
          hpXsc = 39.0+75*(pLab - 1.2)/(std::pow(pLab,3.0) + 0.15);

          fElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
        }
        fTotalXsc = hpXsc;
      }
      if( neutron )     // pn to be np
      {
        if( pLab < 0.02 )
        {
          hnXsc = 4100+30*std::pow(G4Log(1.3/pLab),3.6); // was as pLab < 0.8
	  fElasticXsc = hnXsc;
        }      
        else if( pLab < 0.8 )
        {
          hnXsc = 33+30*std::pow(G4Log(pLab/1.3),4.0);
	  fElasticXsc = hnXsc;
        }      
        else if( pLab < 1.05 )
        {
          hnXsc = 33+30*std::pow(G4Log(pLab/0.95),2.0);
          fElasticXsc =  6 + 52/( G4Log(0.511/pLab)*G4Log(0.511/pLab) + 1.6 );
        }
        else if( pLab < 1.4 )
        {
          hnXsc = 33+30*std::pow(G4Log(pLab/0.95),2.0);
          fElasticXsc =  6 + 52/( G4Log(0.511/pLab)*G4Log(0.511/pLab) + 1.6 );
        }
        else    // 1.4 < pLab < 10.  )
        {
          hnXsc = 33.3 + 20.8*(std::pow(pLab,2.0) - 1.35)/(std::pow(pLab,2.50) + 0.95);
          
          fElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
        }
        fTotalXsc = hnXsc;
      }
    }    
  } 
  else if( theParticle == theAProton && pORn ) /////////////////// p_bar ///////////////////////////
  {
    if( proton )
    {
      xsection  = 35.45 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) + 33.34*std::pow(sMand,-eta2);
    }
    if( neutron ) // ???
    {
      xsection = 35.80 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) + 30.*std::pow(sMand,-eta2);
    }
    fTotalXsc = xsection;
  } 
  else if( theParticle == thePiPlus && pORn ) // pi+ /////////////////////////////////////////////
  {
    if( proton ) // pi+ p
    {
      if( pLab < 0.28 )
      {
        hpXsc       = 10./((logP + 1.273)*(logP + 1.273) + 0.05);
        fElasticXsc = hpXsc;
      }
      else if( pLab < 0.4 )
      {
        hpXsc       = 14./( (logP + 1.273)*(logP + 1.273) + 0.07);
        fElasticXsc = hpXsc;
      }
      else if( pLab < 0.68 )
      {
        hpXsc       = 14./( (logP + 1.273)*(logP + 1.273) + 0.07);
        fElasticXsc = hpXsc;
      }
      else if( pLab < 0.85 )
      {
        G4double Ex4 = 88*(G4Log(pLab/0.77))*(G4Log(pLab/0.77));
        hpXsc        = Ex4 + 14.9;
        fElasticXsc = hpXsc*G4Exp(-3.*(pLab - 0.68));  
      }
      else if( pLab < 1.15 )
      {
        G4double Ex4 = 88*(G4Log(pLab/0.77))*(G4Log(pLab/0.77));
        hpXsc        = Ex4 + 14.9;

        fElasticXsc = 6.0 + 1.4/(( pLab - 1.4)*( pLab - 1.4) + 0.1);
      }
      else if( pLab < 1.4) // ns original
      {
        G4double Ex1 = 3.2*G4Exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
        G4double Ex2 = 12*G4Exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
        hpXsc        = Ex1 + Ex2 + 27.5;
        fElasticXsc = 6.0 + 1.4/(( pLab - 1.4)*( pLab - 1.4) + 0.1);
      }
      else if( pLab < 2.0 ) // ns original
      {
        G4double Ex1 = 3.2*G4Exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
        G4double Ex2 = 12*G4Exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
        hpXsc        = Ex1 + Ex2 + 27.5;
        fElasticXsc = 3.0 + 1.36/( (logP - 0.336)*(logP - 0.336) + 0.08);    
      }
      else if( pLab < 3.5 ) // ns original
      {
        G4double Ex1 = 3.2*G4Exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
        G4double Ex2 = 12*G4Exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
        hpXsc        = Ex1 + Ex2 + 27.5;
        fElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
      }
      else if( pLab < 200. ) // my
      {
        hpXsc = 10.6 + 2.*G4Log(pE) + 25*std::pow(pE, -0.43 ); // ns original
        // hpXsc = GetHadronNucleonXscPDG(aParticle, nucleon )/millibarn;
        fElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
      }
      else //  pLab > 100 // my
      {
        hpXsc = GetHadronNucleonXscPDG(aParticle, nucleon )/millibarn;
        fElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
      }
      fTotalXsc = hpXsc;
    }    
    if( neutron )  // pi+ n = pi- p??
    {
      if( pLab < 0.28 ) 
      {
        hnXsc       = 0.288/((pLab - 0.28)*(pLab - 0.28) + 0.004);
        fElasticXsc = 1.8/((logP + 1.273)*(logP + 1.273) + 0.07);
      }
      else if( pLab < 0.395676 ) // first peak
      {
        hnXsc       = 0.648/((pLab - 0.28)*(pLab - 0.28) + 0.009);
        fElasticXsc = 0.257/((pLab - 0.28)*(pLab - 0.28) + 0.01);
       }
      else if( pLab < 0.5 )
      {
        hnXsc       = 26 + 110*(G4Log(pLab/0.48))*(G4Log(pLab/0.48));
        fElasticXsc = 0.37*hnXsc;
      }
      else if( pLab < 0.65 )
      {
        hnXsc       = 26 + 110*(G4Log(pLab/0.48))*(G4Log(pLab/0.48));
        fElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
      }
      else if( pLab < 0.72 )
      {
        hnXsc = 36.1 + 10*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
      }
      else if( pLab < 0.88 )
      {
        hnXsc = 36.1 + 10*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
      }
      else if( pLab < 1.03 )
      {
        hnXsc = 36.1 + 10*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fElasticXsc = 2.0 + 0.4/((pLab - 1.03)*(pLab - 1.03) + 0.016);
      }
      else if( pLab < 1.15 )
      {
        hnXsc = 36.1 + 10*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fElasticXsc = 2.0 + 0.4/((pLab - 1.03)*(pLab - 1.03) + 0.016);
      }
      else if( pLab < 1.3 )
      {
        hnXsc = 36.1 + 10*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fElasticXsc = 3. + 13./pLab;
      }
      else if( pLab < 2.6 ) // < 3.0) // ns original
      {
        hnXsc = 36.1 + 0.079-4.313*G4Log(pLab)+
                3*G4Exp(-(pLab-2.1)*(pLab-2.1)/0.4/0.4)+
                1.5*G4Exp(-(pLab-1.4)*(pLab-1.4)/0.12/0.12);
        fElasticXsc = 3. + 13./pLab; 
      }
      else if( pLab < 20. ) // < 3.0) // ns original
      {
        hnXsc = 36.1 + 0.079 - 4.313*G4Log(pLab)+
                3*G4Exp(-(pLab-2.1)*(pLab-2.1)/0.4/0.4)+
                1.5*G4Exp(-(pLab-1.4)*(pLab-1.4)/0.12/0.12);
        fElasticXsc = 3. + 13./pLab; 
      }
      else   // mb 
      {
        hnXsc = GetHadronNucleonXscPDG(aParticle, nucleon )/millibarn;
        fElasticXsc = 3. + 13./pLab;
      }
      fTotalXsc = hnXsc;
    }
  } 
  else if( theParticle == thePiMinus && pORn ) /// pi- ////////////////////////////////////////////
  {
    if( neutron )     // pi- n = pi+ p??
    {
      if( pLab < 0.28 )
      {
        hnXsc       = 10./((logP + 1.273)*(logP + 1.273) + 0.05);
        fElasticXsc = hnXsc;
      }
      else if( pLab < 0.4 )
      {
        hnXsc       = 14./( (logP + 1.273)*(logP + 1.273) + 0.07);
        fElasticXsc = hnXsc;
      }
      else if( pLab < 0.68 )
      {
        hnXsc       = 14./( (logP + 1.273)*(logP + 1.273) + 0.07);
        fElasticXsc = hnXsc;
      }
      else if( pLab < 0.85 )
      {
        G4double Ex4 = 88*(G4Log(pLab/0.77))*(G4Log(pLab/0.77));
        hnXsc        = Ex4 + 14.9;
        fElasticXsc = hnXsc*G4Exp(-3.*(pLab - 0.68));  
      }
      else if( pLab < 1.15 )
      {
        G4double Ex4 = 88*(G4Log(pLab/0.77))*(G4Log(pLab/0.77));
        hnXsc        = Ex4 + 14.9;

        fElasticXsc = 6.0 + 1.4/(( pLab - 1.4)*( pLab - 1.4) + 0.1);
      }
      else if( pLab < 1.4) // ns original
      {
        G4double Ex1 = 3.2*G4Exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
        G4double Ex2 = 12*G4Exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
        hnXsc        = Ex1 + Ex2 + 27.5;
        fElasticXsc = 6.0 + 1.4/(( pLab - 1.4)*( pLab - 1.4) + 0.1);
      }
      else if( pLab < 2.0 ) // ns original
      {
        G4double Ex1 = 3.2*G4Exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
        G4double Ex2 = 12*G4Exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
        hnXsc        = Ex1 + Ex2 + 27.5;
        fElasticXsc = 3.0 + 1.36/( (logP - 0.336)*(logP - 0.336) + 0.08);    
      }
      else if( pLab < 3.5 ) // ns original
      {
        G4double Ex1 = 3.2*G4Exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
        G4double Ex2 = 12*G4Exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
        hnXsc        = Ex1 + Ex2 + 27.5;
        fElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
      }
      else if( pLab < 200. ) // my
      {
        hnXsc = 10.6 + 2.*G4Log(pE) + 25*std::pow(pE, -0.43 ); // ns original
        fElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
      }
      else //  pLab > 100 // my
      {
        hnXsc = GetHadronNucleonXscPDG(aParticle, nucleon )/millibarn;
        fElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
      }
      fTotalXsc = hnXsc;
    }
    if( proton )    // pi- p
    {
      if( pLab < 0.28 ) 
      {
        hpXsc       = 0.288/((pLab - 0.28)*(pLab - 0.28) + 0.004);
        fElasticXsc = 1.8/((logP + 1.273)*(logP + 1.273) + 0.07);
      }
      else if( pLab < 0.395676 ) // first peak
      {
        hpXsc       = 0.648/((pLab - 0.28)*(pLab - 0.28) + 0.009);
        fElasticXsc = 0.257/((pLab - 0.28)*(pLab - 0.28) + 0.01);
       }
      else if( pLab < 0.5 )
      {
        hpXsc       = 26 + 110*(G4Log(pLab/0.48))*(G4Log(pLab/0.48));
        fElasticXsc = 0.37*hpXsc;
      }
      else if( pLab < 0.65 )
      {
        hpXsc       = 26 + 110*(G4Log(pLab/0.48))*(G4Log(pLab/0.48));
        fElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
      }
      else if( pLab < 0.72 )
      {
        hpXsc = 36.1+
                10*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
      }
      else if( pLab < 0.88 )
      {
        hpXsc = 36.1+
                10*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
      }
      else if( pLab < 1.03 )
      {
        hpXsc = 36.1+
                10*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fElasticXsc = 2.0 + 0.4/((pLab - 1.03)*(pLab - 1.03) + 0.016);
      }
      else if( pLab < 1.15 )
      {
        hpXsc = 36.1+
                10*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fElasticXsc = 2.0 + 0.4/((pLab - 1.03)*(pLab - 1.03) + 0.016);
      }
      else if( pLab < 1.3 )
      {
        hpXsc = 36.1+
                10*G4Exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*G4Exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fElasticXsc = 3. + 13./pLab;
      }
      else if( pLab < 2.6 ) // < 3.0) // ns original
      {
        hpXsc = 36.1+0.079-4.313*G4Log(pLab)+
                3*G4Exp(-(pLab-2.1)*(pLab-2.1)/0.4/0.4)+
                1.5*G4Exp(-(pLab-1.4)*(pLab-1.4)/0.12/0.12);
        fElasticXsc = 3. +13./pLab; // *G4Log(pLab*6.79);
      }
      else   // mb
      {
        hpXsc = GetHadronNucleonXscPDG(aParticle, nucleon )/millibarn;
        fElasticXsc = 3. + 13./pLab;
      }
      fTotalXsc = hpXsc;
    }
  } 
  else if( (theParticle == theKMinus || theParticle == theK0S) && proton )   // Kmp/K0p /////////////////////////////////
  {
    if( pLab < pMin)
    {
      G4double psp = pLab*std::sqrt(pLab);
      fElasticXsc  = 5.2/psp;
      fTotalXsc    = 14./psp;
    }
    else if( pLab > pMax )
    {
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc           = cofLogE*ld2 + 2.23;
      fTotalXsc           = 1.1*cofLogT*ld2 + 19.7;
    }
    else
    {
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      G4double sp  = std::sqrt(pLab);
      G4double psp = pLab*sp;
      G4double p2  = pLab*pLab;
      G4double p4  = p2*p2;
      G4double lm  = pLab - .39;
      G4double md  = lm*lm + .000356;

      G4double lh1  = pLab - 0.78;
      G4double hd1  = lh1*lh1 + .00166;
 
      G4double lh  = pLab - 1.01;
      G4double hd  = lh*lh + .011;

      G4double lh2  = pLab - 1.63;
      G4double hd2  = lh2*lh2 + .007;

      fElasticXsc  = 5.2/psp + (1.1*cofLogE*ld2 + 2.23)/(1. - .7/sp + .075/p4) 
	             + .004/md + 0.005/hd1+ 0.01/hd2 +.15/hd; // small peaks were added

      fTotalXsc    = 14./psp + (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4) 
	             + .006/md  + 0.01/hd1+ 0.02/hd2 + .20/hd ;
    }
  }
  else if( (theParticle == theKMinus || theParticle == theK0S) && neutron )   // Kmn/K0n //////////////////////////////
  {
    if( pLab > pMax )
    {
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc           = cofLogE*ld2 + 2.23;
      fTotalXsc           = 1.1*cofLogT*ld2 + 19.7;
    }
    else
    {
 
      G4double lh  = pLab - 0.98;
      G4double hd  = lh*lh + .021;

      G4double LogPlab = G4Log( pLab );
      G4double sqrLogPlab = LogPlab * LogPlab;

      fElasticXsc  = // 5.2/psp + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .075/p4) + .004/md 
                     5.0 +  8.1*std::pow(pLab,-1.8 ) + 0.16*sqrLogPlab - 1.3*LogPlab + .15/hd;
      fTotalXsc    = // 14./psp + 
                     //  (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4) 
	// WP                     25.2 +  0. *std::pow(pLab, 0.  ) + 0.38*sqrLogPlab - 2.9*LogPlab	             
                     25.2 +  0.38*sqrLogPlab - 2.9*LogPlab	             
                     //       + .006/md  + 0.01/hd1+ 0.02/hd2 
                        + 0.30/hd ;
    }
  }
  else if(  (theParticle == theKPlus || theParticle == theK0L) && proton  )  // Kpp/aKp ////////////////////////
  {
    if( pLab < pMin )
    {
      G4double lr = pLab - .38;
      G4double lm = pLab - 1.;
      G4double md = lm*lm + .392;   
      fElasticXsc = .7/(lr*lr + .076) + 2./md;
      fTotalXsc   = .7/(lr*lr + .076) + 2.6/md;
    }
    else if( pLab > pMax )
    {
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc           = cofLogE*ld2 + 2.23;
      fTotalXsc           = cofLogT*ld2 + 19.2;
    }
    else
    {
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      G4double lr  = pLab - .38;
      G4double LE  = .7/(lr*lr + .076);
      G4double sp  = std::sqrt(pLab);
      G4double p2  = pLab*pLab;
      G4double p4  = p2*p2;
      G4double lm  = pLab - 1.;
      G4double md  = lm*lm + .392;
      fElasticXsc  = LE + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
      fTotalXsc    = LE + (cofLogT*ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 2.6/md;
    }
  }
  else if(  (theParticle == theKPlus || theParticle == theK0L) && neutron  )  // Kpn/aKn ///////////////////////
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
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc           = cofLogE*ld2 + 2.23;
      fTotalXsc           = cofLogT*ld2 + 19.2;
    }
    else
    {
      G4double ld  = G4Log(pLab) - minLogP;
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
  else if( theParticle == theSMinus && pORn ) 
  {
    xsection  = 35.20 + B*std::pow(G4Log(sMand/s0),2.) 
                          - 199.*std::pow(sMand,-eta1) + 264.*std::pow(sMand,-eta2);
  } 
  else if( theParticle == theGamma && pORn ) // modify later on
  {
    xsection  = 0.0 + B*std::pow(G4Log(sMand/s0),2.) 
      + 0.032*std::pow(sMand,-eta1); // WP - 0.0*std::pow(sMand,-eta2);
    fTotalXsc = xsection;   
  } 
  else  // other then p,n,pi+,pi-,K+,K- as proton ??? 
  {
    if( proton )
    {
      xsection  = 35.45 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2);
    }
    if( neutron )
    {
      xsection += 35.80 + B*std::pow(G4Log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2);
    }
    fTotalXsc = xsection;
  } 
  fTotalXsc   *= millibarn; // parametrised in mb
  fElasticXsc *= millibarn; // parametrised in mb

  if( proton && aParticle->GetDefinition()->GetPDGCharge() > 0. )
  {
    G4double cB = GetCoulombBarrier(aParticle, nucleon);
    fTotalXsc   *= cB;
    fElasticXsc *= cB; 
  }
  fInelasticXsc = fTotalXsc - fElasticXsc;
  if( fInelasticXsc < 0. ) fInelasticXsc = 0.;

  // G4cout<<fTotalXsc/millibarn<<"; "<<fElasticXsc/millibarn<<"; "<<fInelasticXsc/millibarn<<G4endl;

  return fTotalXsc;
}

/////////////////////////////////////////////////////////////////////////////////////
//
// Returns kaon-nucleon cross-section based on smoothed NS for GG model

G4double 
G4HadronNucleonXsc::GetKaonNucleonXscGG(const G4DynamicParticle* aParticle, 
                                          const G4ParticleDefinition* nucleon  )
{
  G4double pLab = aParticle->GetMomentum().mag();

  pLab /= GeV;
  G4double LogPlab = G4Log( pLab );
  G4double sqrLogPlab = LogPlab * LogPlab;

  G4double minLogP = 3.5;       // min of (lnP-minLogP)^2 
  G4double cofLogE = .0557;     // elastic (lnP-minLogP)^2 
  G4double cofLogT = .3;        // total (lnP-minLogP)^2 
  G4double pMin = .1;        // fast LE calculation 
  G4double pMax = 1000.;     // fast HE calculation 

  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();

  G4bool proton = (nucleon == theProton);
  G4bool neutron = (nucleon == theNeutron);

  if(  (theParticle == theKMinus || theParticle == theK0S) && proton ) // (K-,K0)on p ////////////////////////////
  {

    if( pLab < pMin)
    {
      G4double psp = pLab*std::sqrt(pLab);
      fElasticXsc  = 5.2/psp;
      fTotalXsc    = 14./psp;
    }
    else if( pLab > pMax )
    {
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc           = cofLogE*ld2 + 2.23;
      fTotalXsc           = 1.1*cofLogT*ld2 + 19.7;
    }
    else
    {
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      G4double sp  = std::sqrt(pLab);
      G4double psp = pLab*sp;
      G4double p2  = pLab*pLab;
      G4double p4  = p2*p2;
 
      G4double lh  = pLab - 0.98;
      G4double hd  = lh*lh + .045;


      fElasticXsc  = 5.2/psp + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .075/p4) // + .004/md 
               + .15/hd;
      fTotalXsc    = 14./psp + (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4) 
	              //  + .006/md  + 0.01/hd1 + 0.02/hd2 
                     + .60/hd;
    }
  }
  else if( (theParticle == theKMinus || theParticle == theK0S) && neutron )   // Kmn/K0n /////////////////////////////
  {
    if( pLab > pMax )
    {
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc           = cofLogE*ld2 + 2.23;
      fTotalXsc           = 1.1*cofLogT*ld2 + 19.7;
    }
    else
    {
 
      G4double lh  = pLab - 0.98;
      G4double hd  = lh*lh + .045;

      fElasticXsc  = // 5.2/psp + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .075/p4) + .004/md 
                     5.0 +  8.1*std::pow(pLab,-1.8 ) + 0.16*sqrLogPlab - 1.3*LogPlab + .15/hd;
      fTotalXsc    = // 14./psp + 
                     //  (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4) 
	// WP                     25.2 +  0. *std::pow(pLab, 0.  ) + 0.38*sqrLogPlab - 2.9*LogPlab	             
                     25.2 + 0.38*sqrLogPlab - 2.9*LogPlab	             
                     //       + .006/md  + 0.01/hd1+ 0.02/hd2 
                        + 0.60/hd ;
    }
  }
  else if(  (theParticle == theKPlus || theParticle == theK0L) && proton )  // Kpp/aKp //////////////////////
  {
    if( pLab < pMin )
    {
      G4double lr = pLab - .38;
      G4double lm = pLab - 1.;
      G4double md = lm*lm + .392;   
      fElasticXsc = .7/(lr*lr + .076) + 2./md;
      fTotalXsc   = // .7/(lr*lr + .076) + 
                2.6/md;
    }
    else if( pLab > pMax )
    {
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc           = cofLogE*ld2 + 2.23;
      fTotalXsc           = cofLogT*ld2 + 19.2;
    }
    else
    {
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      G4double lr  = pLab - .38;
      G4double LE  = .7/(lr*lr + .076);
      G4double sp  = std::sqrt(pLab);
      G4double p2  = pLab*pLab;
      G4double p4  = p2*p2;
      G4double lm  = pLab - 0.8;
      G4double md  = lm*lm + .652;
      fElasticXsc  = LE + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
      fTotalXsc    = (cofLogT*ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 7.6/md; // + LE;
    }
  }
  else if( (theParticle == theKPlus || theParticle == theK0L) && neutron )  // Kpn/aKn //////////////////////////////////
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
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      fElasticXsc           = cofLogE*ld2 + 2.23;
      fTotalXsc           = cofLogT*ld2 + 19.2;
    }
    else
    {
      G4double ld  = G4Log(pLab) - minLogP;
      G4double ld2 = ld*ld;
      G4double sp  = std::sqrt(pLab);
      G4double p2  = pLab*pLab;
      G4double p4  = p2*p2;
      G4double lm  = pLab - 0.8;
      G4double md  = lm*lm + .652;
      fElasticXsc  = (cofLogE*ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
      fTotalXsc    = (cofLogT*ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 7.6/md;
    }
  }
  fTotalXsc   *= millibarn; // parametrised in mb
  fElasticXsc *= millibarn; // parametrised in mb

  if( proton && aParticle->GetDefinition()->GetPDGCharge() > 0. )
  {
    G4double cB = GetCoulombBarrier(aParticle, nucleon);
    fTotalXsc   *= cB;
    fElasticXsc *= cB; 
  }
  fInelasticXsc = fTotalXsc - fElasticXsc;
  if( fInelasticXsc < 0. ) fInelasticXsc = 0.;

  // G4cout<<fTotalXsc/millibarn<<"; "<<fElasticXsc/millibarn<<"; "<<fInelasticXsc/millibarn<<G4endl;

  return fTotalXsc;
}

/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon cross-section based on V. Uzjinsky parametrisation of
// data from G4FTFCrossSection class

G4double 
G4HadronNucleonXsc::GetHadronNucleonXscVU(const G4DynamicParticle* aParticle, 
                                          const G4ParticleDefinition* nucleon  )
{
  G4int PDGcode = aParticle->GetDefinition()->GetPDGEncoding();
  G4int absPDGcode = std::abs(PDGcode);
  G4double Elab = aParticle->GetTotalEnergy();              
                          // (s - 2*0.88*GeV*GeV)/(2*0.939*GeV)/GeV;
  G4double Plab = aParticle->GetMomentum().mag();            
                          // std::sqrt(Elab * Elab - 0.88);

  Elab /= GeV;
  Plab /= GeV;

  G4double LogPlab = G4Log( Plab );
  G4double sqrLogPlab = LogPlab * LogPlab;

  G4bool pORn = (nucleon == theProton || nucleon == theNeutron  );  
  G4bool proton = (nucleon == theProton);
  G4bool neutron = (nucleon == theNeutron);

   
  if( absPDGcode > 1000 && pORn )  //------Projectile is baryon -
  {
    if(proton)
    {
      // WP      fTotalXsc   = 48.0 +  0. *std::pow(Plab, 0.  ) + 0.522*sqrLogPlab - 4.51*LogPlab;
      fTotalXsc   = 48.0 + 0.522*sqrLogPlab - 4.51*LogPlab;
      fElasticXsc = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;
    }
    if(neutron)
    {    
      // WP      fTotalXsc   = 47.3 +  0. *std::pow(Plab, 0.  ) + 0.513*sqrLogPlab - 4.27*LogPlab;
      fTotalXsc   = 47.3  + 0.513*sqrLogPlab - 4.27*LogPlab;
      fElasticXsc = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;
    }
  }
  else if( PDGcode ==  211  && pORn )  //------Projectile is PionPlus ----
  {
    if(proton)
    {
      fTotalXsc  = 16.4 + 19.3 *std::pow(Plab,-0.42) + 0.19 *sqrLogPlab - 0.0 *LogPlab;
      fElasticXsc =  0.0 + 11.4*std::pow(Plab,-0.40) + 0.079*sqrLogPlab - 0.0 *LogPlab;
    }
    if(neutron)
    {    
      fTotalXsc   =  33.0 + 14.0 *std::pow(Plab,-1.36) + 0.456*sqrLogPlab - 4.03*LogPlab;
      fElasticXsc = 1.76 + 11.2*std::pow(Plab,-0.64) + 0.043*sqrLogPlab - 0.0 *LogPlab;
    }
  }
  else if( PDGcode == -211  && pORn )  //------Projectile is PionMinus ----
  {
    if(proton)
    {
      fTotalXsc   = 33.0 + 14.0 *std::pow(Plab,-1.36) + 0.456*sqrLogPlab - 4.03*LogPlab;
      fElasticXsc = 1.76 + 11.2*std::pow(Plab,-0.64) + 0.043*sqrLogPlab - 0.0 *LogPlab;
    }
    if(neutron)
    {    
      fTotalXsc   = 16.4 + 19.3 *std::pow(Plab,-0.42) + 0.19 *sqrLogPlab - 0.0 *LogPlab;
      fElasticXsc =  0.0 + 11.4*std::pow(Plab,-0.40) + 0.079*sqrLogPlab - 0.0 *LogPlab;
    }
  }
  else if( PDGcode ==  111  && pORn )  //------Projectile is PionZero  --
  {
    if(proton)
    {
      fTotalXsc   = (16.4 + 19.3 *std::pow(Plab,-0.42) + 0.19 *sqrLogPlab - 0.0 *LogPlab +   //Pi+
                        33.0 + 14.0 *std::pow(Plab,-1.36) + 0.456*sqrLogPlab - 4.03*LogPlab)/2; //Pi-

      fElasticXsc = ( 0.0 + 11.4*std::pow(Plab,-0.40) + 0.079*sqrLogPlab - 0.0 *LogPlab +    //Pi+
                         1.76 + 11.2*std::pow(Plab,-0.64) + 0.043*sqrLogPlab - 0.0 *LogPlab)/2; //Pi-

    }
    if(neutron)
    {    
      fTotalXsc   = (33.0 + 14.0 *std::pow(Plab,-1.36) + 0.456*sqrLogPlab - 4.03*LogPlab +   //Pi+
                        16.4 + 19.3 *std::pow(Plab,-0.42) + 0.19 *sqrLogPlab - 0.0 *LogPlab)/2; //Pi-
      fElasticXsc = ( 1.76 + 11.2*std::pow(Plab,-0.64) + 0.043*sqrLogPlab - 0.0 *LogPlab +   //Pi+
                         0.0  + 11.4*std::pow(Plab,-0.40) + 0.079*sqrLogPlab - 0.0 *LogPlab)/2; //Pi-
    }
  }
  else if( PDGcode == 321  && pORn )    //------Projectile is KaonPlus --
  {
    if(proton)
    {
      // WP      fTotalXsc   = 18.1 +  0. *std::pow(Plab, 0.  ) + 0.26 *sqrLogPlab - 1.0 *LogPlab;
      fTotalXsc   = 18.1 +  0.26 *sqrLogPlab - 1.0 *LogPlab;
      fElasticXsc =  5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab - 1.3 *LogPlab;
    }
    if(neutron)
    {    
      // WP      fTotalXsc   = 18.7 +  0. *std::pow(Plab, 0.  ) + 0.21 *sqrLogPlab - 0.89*LogPlab;
      // WP      fElasticXsc =  7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab - 2.4 *LogPlab;
      fTotalXsc   = 18.7  + 0.21 *sqrLogPlab - 0.89*LogPlab;
      fElasticXsc =  7.3  + 0.29 *sqrLogPlab - 2.4 *LogPlab;
    }
  }
  else if( PDGcode ==-321  && pORn )  //------Projectile is KaonMinus ----
  {
    if(proton)
    {
      // WP      fTotalXsc   = 32.1 +  0. *std::pow(Plab, 0.  ) + 0.66*sqrLogPlab - 5.6*LogPlab;
      // WP      fElasticXsc =  7.3 +  0. *std::pow(Plab,-0.  ) + 0.29*sqrLogPlab - 2.4*LogPlab;
      fTotalXsc   = 32.1 + 0.66*sqrLogPlab - 5.6*LogPlab;
      fElasticXsc =  7.3 + 0.29*sqrLogPlab - 2.4*LogPlab;
    }
    if(neutron)
    {    
      // WP      fTotalXsc   = 25.2 +  0. *std::pow(Plab, 0.  ) + 0.38*sqrLogPlab - 2.9*LogPlab;
      fTotalXsc   = 25.2 + 0.38*sqrLogPlab - 2.9*LogPlab;
      fElasticXsc =  5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16*sqrLogPlab - 1.3*LogPlab;
    }
  }
  else if( PDGcode == 311  && pORn )  //------Projectile is KaonZero -----
  {
    if(proton)
    {
      // WP      fTotalXsc   = ( 18.1 +  0. *std::pow(Plab, 0.  ) + 0.26 *sqrLogPlab - 1.0 *LogPlab +   //K+
      // WP                      32.1 +  0. *std::pow(Plab, 0.  ) + 0.66 *sqrLogPlab - 5.6 *LogPlab)/2; //K-
      fTotalXsc   = ( 18.1 + 0.26 *sqrLogPlab - 1.0 *LogPlab +   //K+
                      32.1 + 0.66 *sqrLogPlab - 5.6 *LogPlab)/2; //K-
      fElasticXsc = (  5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab - 1.3 *LogPlab +   //K+
                         7.3 + 0.29 *sqrLogPlab - 2.4 *LogPlab)/2; //K-
      // WP                         7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab - 2.4 *LogPlab)/2; //K-
    }
    if(neutron)
    {    
      // WP      fTotalXsc   = ( 18.7 +  0. *std::pow(Plab, 0.  ) + 0.21 *sqrLogPlab - 0.89*LogPlab +   //K+
      // WP                      25.2 +  0. *std::pow(Plab, 0.  ) + 0.38 *sqrLogPlab - 2.9 *LogPlab)/2; //K-
      fTotalXsc   = ( 18.7 + 0.21 *sqrLogPlab - 0.89*LogPlab +   //K+
                      25.2 + 0.38 *sqrLogPlab - 2.9 *LogPlab)/2; //K-
      // WP      fElasticXsc = (  7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab - 2.4 *LogPlab +   //K+
      fElasticXsc = (  7.3 + 0.29 *sqrLogPlab - 2.4 *LogPlab +   //K+
                       5.0 + 8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab - 1.3 *LogPlab)/2; //K-
    }
  }
  else  //------Projectile is undefined, Nucleon assumed
  {
    if(proton)
    {
      // WP      fTotalXsc   = 48.0 +  0. *std::pow(Plab, 0.  ) + 0.522*sqrLogPlab - 4.51*LogPlab;
      fTotalXsc   = 48.0 + 0.522*sqrLogPlab - 4.51*LogPlab;
      fElasticXsc = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;
    }
    if(neutron)
    {    
      // WP      fTotalXsc   = 47.3 +  0. *std::pow(Plab, 0.  ) + 0.513*sqrLogPlab - 4.27*LogPlab;
      fTotalXsc   = 47.3 + 0.513*sqrLogPlab - 4.27*LogPlab;
      fElasticXsc = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;
    }
  }
  fTotalXsc   *= millibarn;
  fElasticXsc *= millibarn;
  fInelasticXsc   = fTotalXsc - fElasticXsc;
  if (fInelasticXsc < 0.) fInelasticXsc = 0.;

  return fTotalXsc;    
}

////////////////////////////////////////////////////////////////////////////////////
//
//

G4double G4HadronNucleonXsc::CalculateEcmValue( const G4double mp , 
                                                const G4double mt , 
                                                const G4double Plab )
{
  G4double Elab = std::sqrt ( mp * mp + Plab * Plab );
  G4double Ecm  = std::sqrt ( mp * mp + mt * mt + 2 * Elab * mt );
  // G4double Pcm  = Plab * mt / Ecm;
  // G4double KEcm = std::sqrt ( Pcm * Pcm + mp * mp ) - mp;

  return Ecm ; // KEcm;
}


////////////////////////////////////////////////////////////////////////////////////
//
//

G4double G4HadronNucleonXsc::CalcMandelstamS( const G4double mp , 
                                                       const G4double mt , 
                                                       const G4double Plab )
{
  G4double Elab = std::sqrt ( mp * mp + Plab * Plab );
  G4double sMand  = mp*mp + mt*mt + 2*Elab*mt ;

  return sMand;
}


///////////////////////////////////////////////////////////////////////////////
//
//

G4double G4HadronNucleonXsc::GetCoulombBarrier(const G4DynamicParticle* aParticle, 
                                               const G4ParticleDefinition* nucleon )
{
  G4double ratio;

  G4double tR = 0.895*fermi, pR;

  if     ( aParticle->GetDefinition() == theProton ) pR = 0.895*fermi;
  else if( aParticle->GetDefinition() == thePiPlus ) pR = 0.663*fermi;
  else if( aParticle->GetDefinition() == theKPlus )  pR = 0.340*fermi;
  else                                               pR = 0.500*fermi;

  G4double pZ = aParticle->GetDefinition()->GetPDGCharge();
  G4double tZ = nucleon->GetPDGCharge();

  G4double pTkin = aParticle->GetKineticEnergy();
  
  G4double pM    = aParticle->GetDefinition()->GetPDGMass(); 
  G4double tM    = nucleon->GetPDGMass();

  G4double pElab = pTkin + pM;

  G4double totEcm  = std::sqrt(pM*pM + tM*tM + 2.*pElab*tM);

  G4double totTcm  = totEcm - pM -tM;

  G4double bC    = fine_structure_const*hbarc*pZ*tZ;
           bC   /= pR + tR;
           bC   /= 2.;  // 4., 2. parametrisation cof ??? vmg

	   // G4cout<<"pTkin = "<<pTkin/GeV<<"; pPlab = "
	   // <<pPlab/GeV<<"; bC = "<<bC/GeV<<"; pTcm = "<<pTcm/GeV<<G4endl;

  if( totTcm <= bC ) ratio = 0.;
  else               ratio = 1. - bC/totTcm;

  // if(ratio < DBL_MIN) ratio = DBL_MIN;
  if( ratio < 0.) ratio = 0.;

  // G4cout <<"ratio = "<<ratio<<G4endl;
  return ratio;
}



/*

////////////////////////////////////////////////////////////////////////////////////
//
// Initialaise K(p,m)-(p,n) total cross section vectors


void G4HadronNucleonXsc::InitialiseKaonNucleonTotXsc()
{
  G4int i = 0, iMax;
  G4double tmpxsc[106];

  // Kp-proton tot xsc

  iMax = 66;
  fKpProtonTotXscVector = G4LPhysicsFreeVector(iMax, fKpProtonTotTkin[0], fKpProtonTotTkin[iMax-1]);
  fKpProtonTotXscVector.SetSpline(true);

  for( i = 0; i < iMax; i++)
  {
    tmpxsc[i] = 0.;

    if( i == 0 || i == iMax-1 ) tmpxsc[i] = fKpProtonTotXsc[i];
    else                        tmpxsc[i] = 0.5*(fKpProtonTotXsc[i-1]+fKpProtonTotXsc[i+1]);

    fKpProtonTotXscVector.PutValues(size_t(i), fKpProtonTotTkin[i], tmpxsc[i]*millibarn);
  }

  // Kp-neutron tot xsc

  iMax = 75;
  fKpNeutronTotXscVector = G4LPhysicsFreeVector(iMax, fKpNeutronTotTkin[0], fKpNeutronTotTkin[iMax-1]);
  fKpNeutronTotXscVector.SetSpline(true);

  for( i = 0; i < iMax; i++)
  {
    tmpxsc[i] = 0.;
    if( i == 0 || i == iMax-1 ) tmpxsc[i] = fKpNeutronTotXsc[i];
    else                        tmpxsc[i] = 0.5*(fKpNeutronTotXsc[i-1]+fKpNeutronTotXsc[i+1]);

    fKpNeutronTotXscVector.PutValues(size_t(i), fKpNeutronTotTkin[i], tmpxsc[i]*millibarn);
  }

  // Km-proton tot xsc

  iMax = 106;
  fKmProtonTotXscVector = G4LPhysicsFreeVector(iMax, fKmProtonTotTkin[0], fKmProtonTotTkin[iMax-1]);
  fKmProtonTotXscVector.SetSpline(true);

  for( i = 0; i < iMax; i++)
  {
    tmpxsc[i] = 0.;

    if( i == 0 || i == iMax-1 ) tmpxsc[i] = fKmProtonTotXsc[i];
    else                        tmpxsc[i] = 0.5*(fKmProtonTotXsc[i-1]+fKmProtonTotXsc[i+1]);

    fKmProtonTotXscVector.PutValues(size_t(i), fKmProtonTotTkin[i], tmpxsc[i]*millibarn);
  }

  // Km-neutron tot xsc

  iMax = 68;
  fKmNeutronTotXscVector = G4LPhysicsFreeVector(iMax, fKmNeutronTotTkin[0], fKmNeutronTotTkin[iMax-1]);
  fKmNeutronTotXscVector.SetSpline(true);

  for( i = 0; i < iMax; i++)
  {
    tmpxsc[i] = 0.;
    if( i == 0 || i == iMax-1 ) tmpxsc[i] = fKmNeutronTotXsc[i];
    else                        tmpxsc[i] = 0.5*(fKmNeutronTotXsc[i-1]+fKmNeutronTotXsc[i+1]);

    fKmNeutronTotXscVector.PutValues(size_t(i), fKmNeutronTotTkin[i], tmpxsc[i]*millibarn);
  }
}

///////////////////////////////////////////////////////
//
// K-nucleon tot xsc (mb) fit data, G4Log(Tkin(MeV))

const G4double G4HadronNucleonXsc::fKpProtonTotXsc[66] = {
0.000000e+00, 1.592400e-01, 3.184700e-01, 7.961800e-01, 1.433120e+00, 2.070060e+00, 
2.866240e+00, 3.582800e+00, 4.378980e+00, 5.015920e+00, 5.573250e+00, 6.449040e+00, 
7.404460e+00, 8.200640e+00, 8.837580e+00, 9.713380e+00, 1.027070e+01, 1.090764e+01, 
1.130573e+01, 1.170382e+01, 1.242038e+01, 1.281847e+01, 1.321656e+01, 1.337580e+01, 
1.345541e+01, 1.329618e+01, 1.265924e+01, 1.242038e+01, 1.250000e+01, 1.305732e+01, 
1.369427e+01, 1.425159e+01, 1.544586e+01, 1.648089e+01, 1.751592e+01, 1.791401e+01, 
1.791401e+01, 1.775478e+01, 1.751592e+01, 1.735669e+01, 1.719745e+01, 1.711783e+01, 
1.703822e+01, 1.695860e+01, 1.695860e+01, 1.695860e+01, 1.695860e+01, 1.687898e+01, 
1.687898e+01, 1.703822e+01, 1.719745e+01, 1.735669e+01, 1.751592e+01, 1.767516e+01, 
1.783439e+01, 1.799363e+01, 1.815287e+01, 1.839172e+01, 1.855095e+01, 1.871019e+01, 
1.886943e+01, 1.918790e+01, 1.942675e+01, 1.966561e+01, 2.006369e+01, 2.054140e+01 
}; // 66


const G4double G4HadronNucleonXsc::fKpProtonTotTkin[66] = {
2.100000e+00, 2.180770e+00, 2.261540e+00, 2.396150e+00, 2.476920e+00, 2.557690e+00, 
2.557690e+00, 2.584620e+00, 2.638460e+00, 2.665380e+00, 2.719230e+00, 2.746150e+00, 
2.800000e+00, 2.853850e+00, 2.934620e+00, 3.042310e+00, 3.150000e+00, 3.311540e+00, 
3.446150e+00, 3.607690e+00, 3.930770e+00, 4.226920e+00, 4.361540e+00, 4.846150e+00, 
4.980770e+00, 5.088460e+00, 5.465380e+00, 5.653850e+00, 5.950000e+00, 6.084620e+00, 
6.246150e+00, 6.300000e+00, 6.380770e+00, 6.515380e+00, 6.730770e+00, 6.838460e+00, 
7.000000e+00, 7.161540e+00, 7.323080e+00, 7.457690e+00, 7.619230e+00, 7.780770e+00, 
7.915380e+00, 8.130770e+00, 8.265380e+00, 8.453850e+00, 8.642310e+00, 8.803850e+00, 
9.019230e+00, 9.234620e+00, 9.530770e+00, 9.773080e+00, 1.001538e+01, 1.017692e+01, 
1.033846e+01, 1.058077e+01, 1.082308e+01, 1.098462e+01, 1.114615e+01, 1.138846e+01, 
1.160385e+01, 1.173846e+01, 1.192692e+01, 1.216923e+01, 1.238461e+01, 1.257308e+01 
}; // 66

const G4double G4HadronNucleonXsc::fKpNeutronTotXsc[75] = {
3.980900e-01, 3.184700e-01, 3.184700e-01, 3.980900e-01, 3.980900e-01, 3.980900e-01, 
3.980900e-01, 3.980900e-01, 3.980900e-01, 4.777100e-01, 3.980900e-01, 3.980900e-01, 
4.777100e-01, 5.573200e-01, 1.035030e+00, 1.512740e+00, 2.149680e+00, 2.786620e+00, 
3.503180e+00, 4.219750e+00, 5.015920e+00, 5.652870e+00, 6.289810e+00, 7.245220e+00, 
8.121020e+00, 8.837580e+00, 9.633760e+00, 1.042994e+01, 1.114650e+01, 1.194268e+01, 
1.265924e+01, 1.329618e+01, 1.393312e+01, 1.449045e+01, 1.496815e+01, 1.552548e+01, 
1.592357e+01, 1.664013e+01, 1.727707e+01, 1.783439e+01, 1.831210e+01, 1.902866e+01, 
1.902866e+01, 1.878981e+01, 1.847134e+01, 1.831210e+01, 1.807325e+01, 1.791401e+01, 
1.783439e+01, 1.767516e+01, 1.759554e+01, 1.743631e+01, 1.743631e+01, 1.751592e+01, 
1.743631e+01, 1.735669e+01, 1.751592e+01, 1.759554e+01, 1.767516e+01, 1.783439e+01, 
1.783439e+01, 1.791401e+01, 1.815287e+01, 1.823248e+01, 1.847134e+01, 1.878981e+01, 
1.894905e+01, 1.902866e+01, 1.934713e+01, 1.966561e+01, 1.990446e+01, 2.014331e+01, 
2.030255e+01, 2.046178e+01, 2.085987e+01 
}; // 75

const G4double G4HadronNucleonXsc::fKpNeutronTotTkin[75] = {
2.692000e-02, 1.615400e-01, 2.961500e-01, 4.576900e-01, 6.461500e-01, 7.538500e-01, 
8.884600e-01, 1.103850e+00, 1.211540e+00, 1.400000e+00, 1.561540e+00, 1.776920e+00, 
1.992310e+00, 2.126920e+00, 2.342310e+00, 2.423080e+00, 2.557690e+00, 2.692310e+00, 
2.800000e+00, 2.988460e+00, 3.203850e+00, 3.365380e+00, 3.500000e+00, 3.688460e+00, 
3.850000e+00, 4.011540e+00, 4.173080e+00, 4.415380e+00, 4.630770e+00, 4.873080e+00, 
5.061540e+00, 5.276920e+00, 5.492310e+00, 5.707690e+00, 5.896150e+00, 6.030770e+00, 
6.138460e+00, 6.219230e+00, 6.273080e+00, 6.326920e+00, 6.407690e+00, 6.650000e+00, 
6.784620e+00, 7.026920e+00, 7.242310e+00, 7.350000e+00, 7.484620e+00, 7.619230e+00, 
7.807690e+00, 7.915380e+00, 8.050000e+00, 8.211540e+00, 8.453850e+00, 8.588460e+00, 
8.830770e+00, 9.073080e+00, 9.288460e+00, 9.476920e+00, 9.665380e+00, 9.826920e+00, 
1.004231e+01, 1.031154e+01, 1.052692e+01, 1.071538e+01, 1.095769e+01, 1.120000e+01, 
1.138846e+01, 1.155000e+01, 1.176538e+01, 1.190000e+01, 1.214231e+01, 1.222308e+01, 
1.238461e+01, 1.246538e+01, 1.265385e+01 
}; // 75

const G4double G4HadronNucleonXsc::fKmProtonTotXsc[106] = {
1.136585e+02, 9.749129e+01, 9.275262e+01, 8.885017e+01, 8.334146e+01, 7.955401e+01, 
7.504530e+01, 7.153658e+01, 6.858537e+01, 6.674913e+01, 6.525784e+01, 6.448781e+01, 
6.360279e+01, 6.255401e+01, 6.127526e+01, 6.032404e+01, 5.997910e+01, 5.443554e+01, 
5.376307e+01, 5.236934e+01, 5.113937e+01, 5.090941e+01, 4.967944e+01, 4.844948e+01, 
4.705575e+01, 4.638327e+01, 4.571080e+01, 4.475958e+01, 4.397213e+01, 4.257840e+01, 
4.102090e+01, 4.090592e+01, 3.906969e+01, 3.839721e+01, 3.756097e+01, 3.644599e+01, 
3.560976e+01, 3.533101e+01, 3.533101e+01, 3.644599e+01, 3.811847e+01, 3.839721e+01, 
3.979094e+01, 4.090592e+01, 4.257840e+01, 4.341463e+01, 4.425087e+01, 4.564460e+01, 
4.759582e+01, 4.703833e+01, 4.843206e+01, 4.787457e+01, 4.452962e+01, 4.202090e+01, 
4.034843e+01, 3.839721e+01, 3.616725e+01, 3.365854e+01, 3.170732e+01, 3.087108e+01, 
3.170732e+01, 3.254355e+01, 3.310104e+01, 3.254355e+01, 3.142857e+01, 3.059233e+01, 
2.947735e+01, 2.891986e+01, 2.836237e+01, 2.752613e+01, 2.696864e+01, 2.641115e+01, 
2.501742e+01, 2.473868e+01, 2.418118e+01, 2.362369e+01, 2.334495e+01, 2.278746e+01, 
2.250871e+01, 2.222997e+01, 2.167247e+01, 2.139373e+01, 2.139373e+01, 2.139373e+01, 
2.111498e+01, 2.083624e+01, 2.055749e+01, 2.083624e+01, 2.055749e+01, 2.083624e+01, 
2.083624e+01, 2.055749e+01, 2.055749e+01, 2.055749e+01, 2.027875e+01, 2.000000e+01, 
2.055749e+01, 2.027875e+01, 2.083624e+01, 2.083624e+01, 2.055749e+01, 2.083624e+01, 
2.083624e+01, 2.083624e+01, 2.139373e+01, 2.139373e+01
}; // 106

const G4double G4HadronNucleonXsc::fKmProtonTotTkin[106] = {
4.017980e+00, 4.125840e+00, 4.179780e+00, 4.251690e+00, 4.287640e+00, 4.341570e+00, 
4.395510e+00, 4.467420e+00, 4.503370e+00, 4.575280e+00, 4.683150e+00, 4.737080e+00, 
4.773030e+00, 4.826970e+00, 4.880900e+00, 4.916850e+00, 4.952810e+00, 4.988760e+00, 
4.988761e+00, 5.006740e+00, 5.006741e+00, 5.042700e+00, 5.078650e+00, 5.114610e+00, 
5.132580e+00, 5.150560e+00, 5.186520e+00, 5.204490e+00, 5.276400e+00, 5.348310e+00, 
5.366290e+00, 5.384270e+00, 5.456180e+00, 5.564040e+00, 5.600000e+00, 5.671910e+00, 
5.743820e+00, 5.833710e+00, 5.905620e+00, 5.977530e+00, 6.085390e+00, 6.085390e+00, 
6.157300e+00, 6.175280e+00, 6.211240e+00, 6.229210e+00, 6.247190e+00, 6.337080e+00, 
6.391010e+00, 6.516850e+00, 6.462920e+00, 6.498880e+00, 6.570790e+00, 6.606740e+00, 
6.660670e+00, 6.678650e+00, 6.696630e+00, 6.732580e+00, 6.804490e+00, 6.876400e+00, 
6.948310e+00, 7.020220e+00, 7.074160e+00, 7.182020e+00, 7.235960e+00, 7.289890e+00, 
7.397750e+00, 7.523600e+00, 7.631460e+00, 7.757300e+00, 7.901120e+00, 8.062920e+00, 
8.260670e+00, 8.386520e+00, 8.530340e+00, 8.674160e+00, 8.817980e+00, 8.943820e+00, 
9.087640e+00, 9.267420e+00, 9.429210e+00, 9.573030e+00, 9.698880e+00, 9.896630e+00, 
1.002247e+01, 1.016629e+01, 1.031011e+01, 1.048989e+01, 1.063371e+01, 1.077753e+01, 
1.095730e+01, 1.108315e+01, 1.120899e+01, 1.135281e+01, 1.149663e+01, 1.162247e+01, 
1.174831e+01, 1.187416e+01, 1.200000e+01, 1.212584e+01, 1.221573e+01, 1.234157e+01, 
1.239551e+01, 1.250337e+01, 1.261124e+01, 1.273708e+01 
}; // 106

const G4double G4HadronNucleonXsc::fKmNeutronTotXsc[68] = {
2.621810e+01, 2.741123e+01, 2.868413e+01, 2.963889e+01, 3.067343e+01, 3.178759e+01, 
3.282148e+01, 3.417466e+01, 3.536778e+01, 3.552620e+01, 3.544576e+01, 3.496756e+01, 
3.433030e+01, 3.401166e+01, 3.313537e+01, 3.257772e+01, 3.178105e+01, 3.138264e+01, 
3.074553e+01, 2.970952e+01, 2.891301e+01, 2.827542e+01, 2.787700e+01, 2.715978e+01, 
2.660181e+01, 2.612394e+01, 2.564574e+01, 2.516721e+01, 2.421098e+01, 2.365235e+01, 
2.317366e+01, 2.261437e+01, 2.237389e+01, 2.205427e+01, 2.181395e+01, 2.165357e+01, 
2.149335e+01, 2.133297e+01, 2.109232e+01, 2.093128e+01, 2.069030e+01, 2.052992e+01, 
2.028927e+01, 2.012824e+01, 1.996737e+01, 1.996590e+01, 1.988530e+01, 1.964432e+01, 
1.948361e+01, 1.940236e+01, 1.940040e+01, 1.931882e+01, 1.947593e+01, 1.947429e+01, 
1.939320e+01, 1.939157e+01, 1.946922e+01, 1.962715e+01, 1.970481e+01, 1.970301e+01, 
1.993958e+01, 2.009669e+01, 2.025380e+01, 2.033178e+01, 2.049003e+01, 2.064747e+01, 
2.080540e+01, 2.096333e+01 
}; // 68

const G4double G4HadronNucleonXsc::fKmNeutronTotTkin[68] = {
5.708500e+00, 5.809560e+00, 5.896270e+00, 5.954120e+00, 5.997630e+00, 6.041160e+00, 
6.142160e+00, 6.171410e+00, 6.272470e+00, 6.344390e+00, 6.416230e+00, 6.459180e+00, 
6.487690e+00, 6.501940e+00, 6.544740e+00, 6.573280e+00, 6.616110e+00, 6.644710e+00, 
6.658840e+00, 6.744700e+00, 6.773150e+00, 6.830410e+00, 6.859010e+00, 6.916240e+00, 
6.973530e+00, 6.987730e+00, 7.030670e+00, 7.102360e+00, 7.173880e+00, 7.288660e+00, 
7.374720e+00, 7.547000e+00, 7.690650e+00, 7.791150e+00, 7.920420e+00, 8.020980e+00, 
8.107160e+00, 8.207720e+00, 8.365740e+00, 8.523790e+00, 8.710560e+00, 8.811110e+00, 
8.969140e+00, 9.127190e+00, 9.270860e+00, 9.400230e+00, 9.486440e+00, 9.673210e+00, 
9.802510e+00, 9.946220e+00, 1.011870e+01, 1.029116e+01, 1.047808e+01, 1.062181e+01, 
1.075114e+01, 1.089488e+01, 1.106739e+01, 1.118244e+01, 1.135496e+01, 1.151307e+01, 
1.171439e+01, 1.190130e+01, 1.208822e+01, 1.223199e+01, 1.231829e+01, 1.247646e+01, 
1.259150e+01, 1.270655e+01 
}; // 68


*/


//
//
///////////////////////////////////////////////////////////////////////////////////////
