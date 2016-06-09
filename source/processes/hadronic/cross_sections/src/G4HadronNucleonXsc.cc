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

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadTmpUtil.hh"


G4HadronNucleonXsc::G4HadronNucleonXsc() 
: fUpperLimit( 10000 * GeV ),
  fLowerLimit( 0.03 * MeV ),
  fTotalXsc(0.0), fElasticXsc(0.0), fInelasticXsc(0.0), fHadronNucleonXsc(0.0)
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
}


G4HadronNucleonXsc::~G4HadronNucleonXsc()
{}


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

    // xsection = At*( 49.51*std::pow(sMand,-0.097) + 0.314*std::log(sMand)*std::log(sMand) );
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
    // xsection = At*( 55.2*std::pow(sMand,-0.255) + 0.346*std::log(sMand)*std::log(sMand) );
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
      xsection = Zt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
		 + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));// on p
    }
    if ( neutron )
    {
      xsection  = Nt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
		      + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2)); // on n pp for nn
    }
  } 
  else if(theParticle == theProton) 
  {
    if ( proton )
    {      
      xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
    }
  } 
  else if(theParticle == theAProton) 
  {
    if ( proton )
    {      
      xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) + 33.34*std::pow(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) + 30.*std::pow(sMand,-eta2));
    }
  } 
  else if(theParticle == thePiPlus && pORn ) 
  {
    xsection  = At*( 20.86 + B*std::pow(std::log(sMand/s0),2.) 
                          + 19.24*std::pow(sMand,-eta1) - 6.03*std::pow(sMand,-eta2));
  } 
  else if(theParticle == thePiMinus && pORn ) 
  {
    xsection  = At*( 20.86 + B*std::pow(std::log(sMand/s0),2.) 
                          + 19.24*std::pow(sMand,-eta1) + 6.03*std::pow(sMand,-eta2));
  } 
  else if(theParticle == theKPlus) 
  {
    if ( proton )
    {      
      xsection  = Zt*( 17.91 + B*std::pow(std::log(sMand/s0),2.) 
                          + 7.14*std::pow(sMand,-eta1) - 13.45*std::pow(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = Nt*( 17.87 + B*std::pow(std::log(sMand/s0),2.) 
                          + 5.17*std::pow(sMand,-eta1) - 7.23*std::pow(sMand,-eta2));
    }
  } 
  else if(theParticle == theKMinus) 
  {
    if ( proton )
    {      
      xsection  = Zt*( 17.91 + B*std::pow(std::log(sMand/s0),2.) 
                          + 7.14*std::pow(sMand,-eta1) + 13.45*std::pow(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = Nt*( 17.87 + B*std::pow(std::log(sMand/s0),2.) 
                          + 5.17*std::pow(sMand,-eta1) + 7.23*std::pow(sMand,-eta2) );
    }
  }
  else if(theParticle == theSMinus && pORn ) 
  {
    xsection  = At*( 35.20 + B*std::pow(std::log(sMand/s0),2.) 
                          - 199.*std::pow(sMand,-eta1) + 264.*std::pow(sMand,-eta2) );
  } 
  else if(theParticle == theGamma && pORn ) // modify later on
  {
    xsection  = At*( 0.0 + B*std::pow(std::log(sMand/s0),2.) 
                          + 0.032*std::pow(sMand,-eta1) - 0.0*std::pow(sMand,-eta2) );
   
  } 
  else  // as proton ??? 
  {
    if ( proton )
    {      
      xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                       + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2) );
    }
    if ( neutron )
    {
      xsection = Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                      + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
    }
  } 
  xsection *= millibarn; // parametrised in mb

  fTotalXsc     = xsection;
  fInelasticXsc = 0.83*xsection;
  fElasticXsc   = fTotalXsc - fInelasticXsc;
  if (fElasticXsc < 0.)fElasticXsc = 0.;

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
  G4double xsection(0), Delta, A0, B0;
  G4int Zt=1, Nt=1, At=1;
  G4double hpXsc(0);
  G4double hnXsc(0);


  G4double targ_mass = 0.939*GeV;  // ~mean neutron and proton ???

  G4double proj_mass     = aParticle->GetMass();
  G4double proj_energy   = aParticle->GetTotalEnergy(); 
  G4double proj_momentum = aParticle->GetMomentum().mag();

  G4double sMand = CalcMandelstamS ( proj_mass , targ_mass , proj_momentum );

  sMand         /= GeV*GeV;  // in GeV for parametrisation
  proj_momentum /= GeV;
  proj_energy   /= GeV;
  proj_mass     /= GeV;

  // General PDG fit constants

  G4double s0   = 5.38*5.38; // in Gev^2
  G4double eta1 = 0.458;
  G4double eta2 = 0.458;
  G4double B    = 0.308;


  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();

  G4bool pORn = (nucleon == theProton || nucleon == theNeutron  );  
  G4bool proton = (nucleon == theProton);
  G4bool neutron = (nucleon == theNeutron);

  if( theParticle == theNeutron && pORn) 
  {
    if( proj_momentum >= 10.)
    // if( proj_momentum >= 2.)
    {
      Delta = 1.;

      if( proj_energy < 40. ) Delta = 0.916+0.0021*proj_energy;

      if(proj_momentum >= 10.)
      {
        B0 = 7.5;
        A0 = 100. - B0*std::log(3.0e7);

        xsection = A0 + B0*std::log(proj_energy) - 11
                  + 103*std::pow(2*0.93827*proj_energy + proj_mass*proj_mass+
                     0.93827*0.93827,-0.165);        //  mb
      }
      fTotalXsc = xsection;
    }
    else
    {
      // nn to be pp

      if(neutron)
      {
        if( proj_momentum < 0.73 )
        {
          hnXsc = 23 + 50*( std::pow( std::log(0.73/proj_momentum), 3.5 ) );
        }
        else if( proj_momentum < 1.05  )
        {
          hnXsc = 23 + 40*(std::log(proj_momentum/0.73))*
                         (std::log(proj_momentum/0.73));
        }
        else  // if( proj_momentum < 10.  )
        {
          hnXsc = 39.0+
              75*(proj_momentum - 1.2)/(std::pow(proj_momentum,3.0) + 0.15);
        }
        fTotalXsc = hnXsc;
      }
      // pn to be np

      if(proton)
      {
        if( proj_momentum < 0.8 )
        {
          hpXsc = 33+30*std::pow(std::log(proj_momentum/1.3),4.0);
        }      
        else if( proj_momentum < 1.4 )
        {
          hpXsc = 33+30*std::pow(std::log(proj_momentum/0.95),2.0);
        }
        else    // if( proj_momentum < 10.  )
        {
          hpXsc = 33.3+
              20.8*(std::pow(proj_momentum,2.0)-1.35)/
                 (std::pow(proj_momentum,2.50)+0.95);
        }
        fTotalXsc = hpXsc;
      }
      // xsection = hpXsc*Zt + hnXsc*Nt;
    }
  } 
  else if(theParticle == theProton && pORn) 
  {
    if( proj_momentum >= 10.)
    // if( proj_momentum >= 2.)
    {
      Delta = 1.;

      if( proj_energy < 40. ) Delta = 0.916+0.0021*proj_energy;

      if(proj_momentum >= 10.)
      {
        B0 = 7.5;
        A0 = 100. - B0*std::log(3.0e7);

        xsection = A0 + B0*std::log(proj_energy) - 11
                  + 103*std::pow(2*0.93827*proj_energy + proj_mass*proj_mass+
                     0.93827*0.93827,-0.165);        //  mb
      }
      fTotalXsc = xsection;
    }
    else
    {
      // pp

      if(proton)
      {
        if( proj_momentum < 0.73 )
        {
          hpXsc = 23 + 50*( std::pow( std::log(0.73/proj_momentum), 3.5 ) );
        }
        else if( proj_momentum < 1.05  )
        {
          hpXsc = 23 + 40*(std::log(proj_momentum/0.73))*
                         (std::log(proj_momentum/0.73));
        }
        else    // if( proj_momentum < 10.  )
        {
           hpXsc = 39.0+
              75*(proj_momentum - 1.2)/(std::pow(proj_momentum,3.0) + 0.15);
        }
        fTotalXsc = hpXsc;
      }
      // pn to be np

      if(neutron)
      {
        if( proj_momentum < 0.8 )
        {
          hnXsc = 33+30*std::pow(std::log(proj_momentum/1.3),4.0);
        }      
        else if( proj_momentum < 1.4 )
        {
          hnXsc = 33+30*std::pow(std::log(proj_momentum/0.95),2.0);
        }
        else   // if( proj_momentum < 10.  )
        {
          hnXsc = 33.3+
              20.8*(std::pow(proj_momentum,2.0)-1.35)/
                 (std::pow(proj_momentum,2.50)+0.95);
        }
        fTotalXsc = hnXsc;
      }
      // xsection = hpXsc*Zt + hnXsc*Nt;
      // xsection = hpXsc*(Zt + Nt);
      // xsection = hnXsc*(Zt + Nt);
    }    
    // xsection *= 0.95;
  } 
  else if(theParticle == theAProton && pORn) 
  {
    if(proton)
    {
      xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) + 33.34*std::pow(sMand,-eta2));
    }
    if(proton)
    {
      xsection = Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) + 30.*std::pow(sMand,-eta2));
    }
    fTotalXsc = xsection;
  } 
  else if(theParticle == thePiPlus && pORn) 
  {
    if(proton)
    {
      if(proj_momentum < 0.4)
      {
        G4double Ex3 = 180*std::exp(-(proj_momentum-0.29)*(proj_momentum-0.29)/0.085/0.085);
        hpXsc      = Ex3+20.0;
      }
      else if(proj_momentum < 1.15)
      {
        G4double Ex4 = 88*(std::log(proj_momentum/0.75))*(std::log(proj_momentum/0.75));
        hpXsc = Ex4+14.0;
      }
      else if(proj_momentum < 3.5)
      {
        G4double Ex1 = 3.2*std::exp(-(proj_momentum-2.55)*(proj_momentum-2.55)/0.55/0.55);
        G4double Ex2 = 12*std::exp(-(proj_momentum-1.47)*(proj_momentum-1.47)/0.225/0.225);
        hpXsc = Ex1+Ex2+27.5;
      }
      else //  if(proj_momentum > 3.5) // mb
      {
        hpXsc = 10.6+2.*std::log(proj_energy)+25*std::pow(proj_energy,-0.43);
      }
      fTotalXsc = hpXsc;
    }    

// pi+n = pi-p??

    if(neutron)
    {
      if(proj_momentum < 0.37)
      {
        hnXsc = 28.0 + 40*std::exp(-(proj_momentum-0.29)*(proj_momentum-0.29)/0.07/0.07);
      }
      else if(proj_momentum<0.65)
      {
        hnXsc = 26+110*(std::log(proj_momentum/0.48))*(std::log(proj_momentum/0.48));
      }
      else if(proj_momentum<1.3)
      {
        hnXsc = 36.1+
                10*std::exp(-(proj_momentum-0.72)*(proj_momentum-0.72)/0.06/0.06)+
                24*std::exp(-(proj_momentum-1.015)*(proj_momentum-1.015)/0.075/0.075);
      }
      else if(proj_momentum<3.0)
      {
        hnXsc = 36.1+0.079-4.313*std::log(proj_momentum)+
                3*std::exp(-(proj_momentum-2.1)*(proj_momentum-2.1)/0.4/0.4)+
                1.5*std::exp(-(proj_momentum-1.4)*(proj_momentum-1.4)/0.12/0.12);
      }
      else   // mb
      {
        hnXsc = 10.6+2*std::log(proj_energy)+30*std::pow(proj_energy,-0.43); 
      }
      fTotalXsc = hnXsc;
    }
    // xsection = hpXsc*Zt + hnXsc*Nt;
  } 
  else if(theParticle == thePiMinus && pORn) 
  {
    // pi-n = pi+p??

    if(neutron)
    {
      if(proj_momentum < 0.4)
      {
        G4double Ex3 = 180*std::exp(-(proj_momentum-0.29)*(proj_momentum-0.29)/0.085/0.085);
        hnXsc      = Ex3+20.0;
      }
      else if(proj_momentum < 1.15)
      {
        G4double Ex4 = 88*(std::log(proj_momentum/0.75))*(std::log(proj_momentum/0.75));
        hnXsc = Ex4+14.0;
      }
      else if(proj_momentum < 3.5)
      {
        G4double Ex1 = 3.2*std::exp(-(proj_momentum-2.55)*(proj_momentum-2.55)/0.55/0.55);
        G4double Ex2 = 12*std::exp(-(proj_momentum-1.47)*(proj_momentum-1.47)/0.225/0.225);
        hnXsc = Ex1+Ex2+27.5;
      }
      else //  if(proj_momentum > 3.5) // mb
      {
      hnXsc = 10.6+2.*std::log(proj_energy)+25*std::pow(proj_energy,-0.43);
      }
      fTotalXsc = hnXsc;
    }
    // pi-p

    if(proton)
    {
      if(proj_momentum < 0.37)
      {
        hpXsc = 28.0 + 40*std::exp(-(proj_momentum-0.29)*(proj_momentum-0.29)/0.07/0.07);
      }
      else if(proj_momentum<0.65)
      {
        hpXsc = 26+110*(std::log(proj_momentum/0.48))*(std::log(proj_momentum/0.48));
      }
      else if(proj_momentum<1.3)
      {
        hpXsc = 36.1+
                10*std::exp(-(proj_momentum-0.72)*(proj_momentum-0.72)/0.06/0.06)+
                24*std::exp(-(proj_momentum-1.015)*(proj_momentum-1.015)/0.075/0.075);
      }
      else if(proj_momentum<3.0)
      {
        hpXsc = 36.1+0.079-4.313*std::log(proj_momentum)+
                3*std::exp(-(proj_momentum-2.1)*(proj_momentum-2.1)/0.4/0.4)+
                1.5*std::exp(-(proj_momentum-1.4)*(proj_momentum-1.4)/0.12/0.12);
      }
      else   // mb
      {
        hpXsc = 10.6+2*std::log(proj_energy)+30*std::pow(proj_energy,-0.43); 
      }
      fTotalXsc = hpXsc;
    }
    // xsection = hpXsc*Zt + hnXsc*Nt;
  } 
  else if(theParticle == theKPlus && pORn) 
  {
    if(proton)
    {
      xsection  = Zt*( 17.91 + B*std::pow(std::log(sMand/s0),2.) 
                          + 7.14*std::pow(sMand,-eta1) - 13.45*std::pow(sMand,-eta2));
    }
    if(neutron)
    {
      xsection = Nt*( 17.87 + B*std::pow(std::log(sMand/s0),2.) 
                          + 5.17*std::pow(sMand,-eta1) - 7.23*std::pow(sMand,-eta2));
    }
    fTotalXsc = xsection;
  } 
  else if(theParticle == theKMinus && pORn) 
  {
    if(proton)
    {
      xsection  = Zt*( 17.91 + B*std::pow(std::log(sMand/s0),2.) 
                          + 7.14*std::pow(sMand,-eta1) + 13.45*std::pow(sMand,-eta2));
    }
    if(neutron)
    {
      xsection = Nt*( 17.87 + B*std::pow(std::log(sMand/s0),2.) 
                          + 5.17*std::pow(sMand,-eta1) + 7.23*std::pow(sMand,-eta2));
    }
    fTotalXsc = xsection;
  }
  else if(theParticle == theSMinus && pORn) 
  {
    xsection  = At*( 35.20 + B*std::pow(std::log(sMand/s0),2.) 
                          - 199.*std::pow(sMand,-eta1) + 264.*std::pow(sMand,-eta2));
  } 
  else if(theParticle == theGamma && pORn) // modify later on
  {
    xsection  = At*( 0.0 + B*std::pow(std::log(sMand/s0),2.) 
                          + 0.032*std::pow(sMand,-eta1) - 0.0*std::pow(sMand,-eta2));
    fTotalXsc = xsection;   
  } 
  else  // as proton ??? 
  {
    if(proton)
    {
      xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2));
    }
    if(neutron)
    {
      xsection += Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
    }
    fTotalXsc = xsection;
  } 
  fTotalXsc *= millibarn; // parametrised in mb
  // xsection  *= millibarn; // parametrised in mb

  fInelasticXsc = 0.83*fTotalXsc;
  fElasticXsc   = fTotalXsc - fInelasticXsc;
  if (fElasticXsc < 0.)fElasticXsc = 0.;

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

  G4double LogPlab = std::log( Plab );
  G4double sqrLogPlab = LogPlab * LogPlab;

  G4bool pORn = (nucleon == theProton || nucleon == theNeutron  );  
  G4bool proton = (nucleon == theProton);
  G4bool neutron = (nucleon == theNeutron);

   
  if( absPDGcode > 1000 && pORn )  //------Projectile is baryon -
  {
    if(proton)
    {
      fTotalXsc   = 48.0 +  0. *std::pow(Plab, 0.  ) + 0.522*sqrLogPlab - 4.51*LogPlab;
      fElasticXsc = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;
    }
    if(neutron)
    {    
      fTotalXsc   = 47.3 +  0. *std::pow(Plab, 0.  ) + 0.513*sqrLogPlab - 4.27*LogPlab;
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
      fTotalXsc   = 18.1 +  0. *std::pow(Plab, 0.  ) + 0.26 *sqrLogPlab - 1.0 *LogPlab;
      fElasticXsc =  5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab - 1.3 *LogPlab;
    }
    if(neutron)
    {    
      fTotalXsc   = 18.7 +  0. *std::pow(Plab, 0.  ) + 0.21 *sqrLogPlab - 0.89*LogPlab;
      fElasticXsc =  7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab - 2.4 *LogPlab;
    }
  }
  else if( PDGcode ==-321  && pORn )  //------Projectile is KaonMinus ----
  {
    if(proton)
    {
      fTotalXsc   = 32.1 +  0. *std::pow(Plab, 0.  ) + 0.66*sqrLogPlab - 5.6*LogPlab;
      fElasticXsc =  7.3 +  0. *std::pow(Plab,-0.  ) + 0.29*sqrLogPlab - 2.4*LogPlab;
    }
    if(neutron)
    {    
      fTotalXsc   = 25.2 +  0. *std::pow(Plab, 0.  ) + 0.38*sqrLogPlab - 2.9*LogPlab;
      fElasticXsc =  5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16*sqrLogPlab - 1.3*LogPlab;
    }
  }
  else if( PDGcode == 311  && pORn )  //------Projectile is KaonZero -----
  {
    if(proton)
    {
      fTotalXsc   = ( 18.1 +  0. *std::pow(Plab, 0.  ) + 0.26 *sqrLogPlab - 1.0 *LogPlab +   //K+
                        32.1 +  0. *std::pow(Plab, 0.  ) + 0.66 *sqrLogPlab - 5.6 *LogPlab)/2; //K-
      fElasticXsc = (  5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab - 1.3 *LogPlab +   //K+
                         7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab - 2.4 *LogPlab)/2; //K-
    }
    if(neutron)
    {    
      fTotalXsc   = ( 18.7 +  0. *std::pow(Plab, 0.  ) + 0.21 *sqrLogPlab - 0.89*LogPlab +   //K+
                         25.2 +  0. *std::pow(Plab, 0.  ) + 0.38 *sqrLogPlab - 2.9 *LogPlab)/2; //K-
      fElasticXsc = (  7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab - 2.4 *LogPlab +   //K+
                         5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab - 1.3 *LogPlab)/2; //K-
    }
  }
  else  //------Projectile is undefined, Nucleon assumed
  {
    if(proton)
    {
      fTotalXsc   = 48.0 +  0. *std::pow(Plab, 0.  ) + 0.522*sqrLogPlab - 4.51*LogPlab;
      fElasticXsc = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;
    }
    if(neutron)
    {    
      fTotalXsc   = 47.3 +  0. *std::pow(Plab, 0.  ) + 0.513*sqrLogPlab - 4.27*LogPlab;
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


//
//
///////////////////////////////////////////////////////////////////////////////////////
