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
// 17.07.06 V. Grichine - first implementation



#include "G4GlauberGribovCrossSection.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"

//////////////////////////////////////////////////////////////////////////////////////
//
//


G4GlauberGribovCrossSection::G4GlauberGribovCrossSection() 
: fUpperLimit( 10000 * GeV ),
  fLowerLimit( 3 * GeV ),
  fRadiusConst( 1.08*fermi )  // 1.1, 1.3 ?
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

///////////////////////////////////////////////////////////////////////////////////////
//
//

G4GlauberGribovCrossSection::~G4GlauberGribovCrossSection()
{
}


////////////////////////////////////////////////////////////////////////////////////////
//
//


G4bool 
G4GlauberGribovCrossSection::IsApplicable(const G4DynamicParticle* aDP, const G4Element*  anElement)
{
  G4bool applicable      = false;
  // G4int baryonNumber     = aDP->GetDefinition()->GetBaryonNumber();
  G4double kineticEnergy = aDP->GetKineticEnergy();

  const G4ParticleDefinition* theParticle = aDP->GetDefinition();
 
  if ( ( kineticEnergy  >= fLowerLimit &&
         anElement->GetZ() > 1. &&      // >=  He
       ( theParticle == theAProton   ||
         theParticle == theGamma     ||
         theParticle == thePiPlus    ||
         theParticle == thePiMinus   ||
         theParticle == theKPlus     ||
         theParticle == theKMinus    || 
         theParticle == theSMinus)      )    ||  

       ( kineticEnergy  >= 0.1*fLowerLimit &&
         anElement->GetZ() > 1. &&      // >=  He
       ( theParticle == theProton    ||
         theParticle == theNeutron      ) )    ) applicable = true;
  return applicable;
}


////////////////////////////////////////////////////////////////////////////////////////
//
// Calculates total and inelastic Xsc, derives elastic as total - inelastic accordong to
// Glauber model with Gribov correction calculated in the dipole approximation on
// light cone. Gaussian density helps to calculate rest integrals of the model.
// [1] B.Z. Kopeliovich, nucl-th/0306044 



G4double G4GlauberGribovCrossSection::
GetCrossSection(const G4DynamicParticle* aParticle, const G4Element* anElement, G4double )
{
  G4double xsection, sigma, cofInelastic, cofTotal, nucleusSquare, ratio;
  G4double R             = GetNucleusRadius(aParticle, anElement); 



  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();

  if(theParticle == theProton || theParticle == theNeutron)
  {
    sigma        = GetHadronNucleaonXscNS(aParticle, anElement);
    cofInelastic = 2.4;
    cofTotal     = 2.0;
  }
  else
  {
    sigma        = GetHadronNucleaonXscPDG(aParticle, anElement);
    cofInelastic = 2.2;
    cofTotal     = 2.0;
  }
  nucleusSquare = cofTotal*pi*R*R;   // basically 2piRR
  ratio = sigma/nucleusSquare;

  xsection =  nucleusSquare*std::log( 1. + ratio );

  fTotalXsc = xsection;

  /*   
  fElasticXsc = 0.5*( xsection - nucleusSquare*ratio/(1.+ratio) );

  if (fElasticXsc < 0.) fElasticXsc = 0.;

  fInelasticXsc = fTotalXsc - fElasticXsc;

  if (fInelasticXsc < 0.) fInelasticXsc = 0.;
  */

  fInelasticXsc = nucleusSquare*std::log( 1. + cofInelastic*ratio )/cofInelastic;
  fElasticXsc   = fTotalXsc - fInelasticXsc;
  if (fElasticXsc < 0.) fElasticXsc = 0.;

  return xsection; 
}

/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to differnt parametrisations:
// [2] E. Levin, hep-ph/9710546
// [3] U. Dersch, et al, hep-ex/9910052
// [4] M.J. Longo, et al, Phys.Rev.Lett. 33 (1974) 725 

G4double 
G4GlauberGribovCrossSection::GetHadronNucleaonXsc(const G4DynamicParticle* aParticle, 
                                                  const G4Element* anElement          )
{
  G4double At = anElement->GetN();  // number of nucleons 
  G4double Zt = anElement->GetZ();  // number of protons


  return GetHadronNucleaonXsc( aParticle, At, Zt );
}




/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to differnt parametrisations:
// [2] E. Levin, hep-ph/9710546
// [3] U. Dersch, et al, hep-ex/9910052
// [4] M.J. Longo, et al, Phys.Rev.Lett. 33 (1974) 725 

G4double 
G4GlauberGribovCrossSection::GetHadronNucleaonXsc(const G4DynamicParticle* aParticle, 
                                                   G4double At,  G4double Zt       )
{
  G4double xsection;


  G4double targ_mass = G4ParticleTable::GetParticleTable()->
  GetIonTable()->GetIonMass( G4int(Zt+0.5) , G4int(At+0.5) );

  targ_mass = 0.939*GeV;  // ~mean neutron and proton ???

  G4double proj_mass     = aParticle->GetMass();
  G4double proj_momentum = aParticle->GetMomentum().mag();
  G4double sMand = CalcMandelstamS ( proj_mass , targ_mass , proj_momentum );

  sMand /= GeV*GeV;  // in GeV for parametrisation
  proj_momentum /= GeV;

  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  

  if(theParticle == theGamma) 
  {
    xsection = At*(0.0677*std::pow(sMand,0.0808) + 0.129*std::pow(sMand,-0.4525));
  } 
  else if(theParticle == theNeutron) // as proton ??? 
  {
    xsection = At*(21.70*std::pow(sMand,0.0808) + 56.08*std::pow(sMand,-0.4525));
  } 
  else if(theParticle == theProton) 
  {
    xsection = At*(21.70*std::pow(sMand,0.0808) + 56.08*std::pow(sMand,-0.4525));
    // xsection = At*( 49.51*std::pow(sMand,-0.097) + 0.314*std::log(sMand)*std::log(sMand) );
    // xsection = At*( 38.4 + 0.85*std::abs(std::pow(log(sMand),1.47)) );
  } 
  else if(theParticle == theAProton) 
  {
    xsection = At*( 21.70*std::pow(sMand,0.0808) + 98.39*std::pow(sMand,-0.4525));
  } 
  else if(theParticle == thePiPlus) 
  {
    xsection = At*(13.63*std::pow(sMand,0.0808) + 27.56*std::pow(sMand,-0.4525));
  } 
  else if(theParticle == thePiMinus) 
  {
    // xsection = At*( 55.2*std::pow(sMand,-0.255) + 0.346*std::log(sMand)*std::log(sMand) );
    xsection = At*(13.63*std::pow(sMand,0.0808) + 36.02*std::pow(sMand,-0.4525));
  } 
  else if(theParticle == theKPlus) 
  {
    xsection = At*(11.82*std::pow(sMand,0.0808) + 8.15*std::pow(sMand,-0.4525));
  } 
  else if(theParticle == theKMinus) 
  {
    xsection = At*(11.82*std::pow(sMand,0.0808) + 26.36*std::pow(sMand,-0.4525));
  }
  else  // as proton ??? 
  {
    xsection = At*(21.70*std::pow(sMand,0.0808) + 56.08*std::pow(sMand,-0.4525));
  } 
  xsection *= millibarn;
  return xsection;
}


/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to PDG parametrisation (2005):
// http://pdg.lbl.gov/2006/reviews/hadronicrpp.pdf

G4double 
G4GlauberGribovCrossSection::GetHadronNucleaonXscPDG(const G4DynamicParticle* aParticle, 
                                                  const G4Element* anElement          )
{
  G4double At = anElement->GetN();  // number of nucleons 
  G4double Zt = anElement->GetZ();  // number of protons


  return GetHadronNucleaonXscPDG( aParticle, At, Zt );
}




/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to PDG parametrisation (2005):
// http://pdg.lbl.gov/2006/reviews/hadronicrpp.pdf
//  At = number of nucleons,  Zt = number of protons 

G4double 
G4GlauberGribovCrossSection::GetHadronNucleaonXscPDG(const G4DynamicParticle* aParticle, 
                                                     G4double At,  G4double Zt )
{
  G4double xsection;

  G4double Nt = At-Zt;              // number of neutrons
  if (Nt < 0.) Nt = 0.;  


  G4double targ_mass = G4ParticleTable::GetParticleTable()->
  GetIonTable()->GetIonMass( G4int(Zt+0.5) , G4int(At+0.5) );

  targ_mass = 0.939*GeV;  // ~mean neutron and proton ???

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
  

  if(theParticle == theNeutron) // proton-neutron fit 
  {
    xsection = Zt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
    xsection  += Nt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
		      + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2)); // pp for nn
  } 
  else if(theParticle == theProton) 
  {
      
      xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2));

      xsection += Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
  } 
  else if(theParticle == theAProton) 
  {
    xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) + 33.34*std::pow(sMand,-eta2));

    xsection += Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) + 30.*std::pow(sMand,-eta2));
  } 
  else if(theParticle == thePiPlus) 
  {
    xsection  = At*( 20.86 + B*std::pow(std::log(sMand/s0),2.) 
                          + 19.24*std::pow(sMand,-eta1) - 6.03*std::pow(sMand,-eta2));
  } 
  else if(theParticle == thePiMinus) 
  {
    xsection  = At*( 20.86 + B*std::pow(std::log(sMand/s0),2.) 
                          + 19.24*std::pow(sMand,-eta1) + 6.03*std::pow(sMand,-eta2));
  } 
  else if(theParticle == theKPlus) 
  {
    xsection  = Zt*( 17.91 + B*std::pow(std::log(sMand/s0),2.) 
                          + 7.14*std::pow(sMand,-eta1) - 13.45*std::pow(sMand,-eta2));

    xsection += Nt*( 17.87 + B*std::pow(std::log(sMand/s0),2.) 
                          + 5.17*std::pow(sMand,-eta1) - 7.23*std::pow(sMand,-eta2));
  } 
  else if(theParticle == theKMinus) 
  {
    xsection  = Zt*( 17.91 + B*std::pow(std::log(sMand/s0),2.) 
                          + 7.14*std::pow(sMand,-eta1) + 13.45*std::pow(sMand,-eta2));

    xsection += Nt*( 17.87 + B*std::pow(std::log(sMand/s0),2.) 
                          + 5.17*std::pow(sMand,-eta1) + 7.23*std::pow(sMand,-eta2));
  }
  else if(theParticle == theSMinus) 
  {
    xsection  = At*( 35.20 + B*std::pow(std::log(sMand/s0),2.) 
                          - 199.*std::pow(sMand,-eta1) + 264.*std::pow(sMand,-eta2));
  } 
  else if(theParticle == theGamma) // modify later on
  {
    xsection  = At*( 0.0 + B*std::pow(std::log(sMand/s0),2.) 
                          + 0.032*std::pow(sMand,-eta1) - 0.0*std::pow(sMand,-eta2));
   
  } 
  else  // as proton ??? 
  {
    xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2));

    xsection += Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
  } 
  xsection *= millibarn; // parametrised in mb
  return xsection;
}
/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon cross-section based on N. Starkov parametrisation of
// data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database

G4double 
G4GlauberGribovCrossSection::GetHadronNucleaonXscNS(const G4DynamicParticle* aParticle, 
                                                  const G4Element* anElement          )
{
  G4double At = anElement->GetN();  // number of nucleons 
  G4double Zt = anElement->GetZ();  // number of protons


  return GetHadronNucleaonXscNS( aParticle, At, Zt );
}




/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon cross-section based on N. Starkov parametrisation of
// data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database

G4double 
G4GlauberGribovCrossSection::GetHadronNucleaonXscNS(const G4DynamicParticle* aParticle, 
                                                     G4double At,  G4double Zt )
{
  G4double xsection, Delta, A0, B0;
  G4double hpXsc(0);
  G4double hnXsc(0);

  G4double Nt = At-Zt;              // number of neutrons
  if (Nt < 0.) Nt = 0.;  


  G4double targ_mass = G4ParticleTable::GetParticleTable()->
  GetIonTable()->GetIonMass( G4int(Zt+0.5) , G4int(At+0.5) );

  targ_mass = 0.939*GeV;  // ~mean neutron and proton ???

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
  

  if(theParticle == theNeutron) 
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
      xsection *= Zt + Nt;
    }
    else
    {
      // nn to be pp

      if( proj_momentum < 10.  )
      {
         hnXsc = 39.0+
              75*(proj_momentum - 1.2)/(std::pow(proj_momentum,3.0) + 0.15);
      }
      if( proj_momentum < 1.05  )
      {
       hnXsc = 23 + 40*(std::log(proj_momentum/0.73))*
                         (std::log(proj_momentum/0.73));
      }
      if( proj_momentum < 0.73 )
      {
        hnXsc = 23 + 50*( std::pow( std::log(0.73/proj_momentum), 3.5 ) );
      }
      // pn to be np

      if( proj_momentum < 10.  )
      {
        hpXsc = 33.3+
              20.8*(std::pow(proj_momentum,2.0)-1.35)/
                 (std::pow(proj_momentum,2.50)+0.95);
      }

      if( proj_momentum < 1.4 )
      {
        hpXsc = 33+30*std::pow(std::log(proj_momentum/0.95),2.0);
      }

      if( proj_momentum < 0.8 )
      {
        hpXsc = 33+30*std::pow(std::log(proj_momentum/1.3),4.0);
      }      
      xsection = hpXsc*Zt + hnXsc*Nt;
    }
  } 
  else if(theParticle == theProton) 
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
      xsection *= Zt + Nt;
    }
    else
    {
      // pp

      if( proj_momentum < 10.  )
      {
         hpXsc = 39.0+
              75*(proj_momentum - 1.2)/(std::pow(proj_momentum,3.0) + 0.15);
      }
      if( proj_momentum < 1.05  )
      {
       hpXsc = 23 + 40*(std::log(proj_momentum/0.73))*
                         (std::log(proj_momentum/0.73));
      }
      if( proj_momentum < 0.73 )
      {
        hpXsc = 23 + 50*( std::pow( std::log(0.73/proj_momentum), 3.5 ) );
      }
      // pn to be np

      if( proj_momentum < 10.  )
      {
        hnXsc = 33.3+
              20.8*(std::pow(proj_momentum,2.0)-1.35)/
                 (std::pow(proj_momentum,2.50)+0.95);
      }

      if( proj_momentum < 1.4 )
      {
        hnXsc = 33+30*std::pow(std::log(proj_momentum/0.95),2.0);
      }

      if( proj_momentum < 0.8 )
      {
        hnXsc = 33+30*std::pow(std::log(proj_momentum/1.3),4.0);
      }      
      xsection = hpXsc*Zt + hnXsc*Nt;
      // xsection = hpXsc*(Zt + Nt);
      // xsection = hnXsc*(Zt + Nt);
    }    
    // xsection *= 0.95;
  } 
  else if(theParticle == theAProton) 
  {
    xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) + 33.34*std::pow(sMand,-eta2));

    xsection += Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) + 30.*std::pow(sMand,-eta2));
  } 
  else if(theParticle == thePiPlus) 
  {
    xsection  = At*( 20.86 + B*std::pow(std::log(sMand/s0),2.) 
                          + 19.24*std::pow(sMand,-eta1) - 6.03*std::pow(sMand,-eta2));
  } 
  else if(theParticle == thePiMinus) 
  {
    xsection  = At*( 20.86 + B*std::pow(std::log(sMand/s0),2.) 
                          + 19.24*std::pow(sMand,-eta1) + 6.03*std::pow(sMand,-eta2));
  } 
  else if(theParticle == theKPlus) 
  {
    xsection  = Zt*( 17.91 + B*std::pow(std::log(sMand/s0),2.) 
                          + 7.14*std::pow(sMand,-eta1) - 13.45*std::pow(sMand,-eta2));

    xsection += Nt*( 17.87 + B*std::pow(std::log(sMand/s0),2.) 
                          + 5.17*std::pow(sMand,-eta1) - 7.23*std::pow(sMand,-eta2));
  } 
  else if(theParticle == theKMinus) 
  {
    xsection  = Zt*( 17.91 + B*std::pow(std::log(sMand/s0),2.) 
                          + 7.14*std::pow(sMand,-eta1) + 13.45*std::pow(sMand,-eta2));

    xsection += Nt*( 17.87 + B*std::pow(std::log(sMand/s0),2.) 
                          + 5.17*std::pow(sMand,-eta1) + 7.23*std::pow(sMand,-eta2));
  }
  else if(theParticle == theSMinus) 
  {
    xsection  = At*( 35.20 + B*std::pow(std::log(sMand/s0),2.) 
                          - 199.*std::pow(sMand,-eta1) + 264.*std::pow(sMand,-eta2));
  } 
  else if(theParticle == theGamma) // modify later on
  {
    xsection  = At*( 0.0 + B*std::pow(std::log(sMand/s0),2.) 
                          + 0.032*std::pow(sMand,-eta1) - 0.0*std::pow(sMand,-eta2));
   
  } 
  else  // as proton ??? 
  {
    xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2));

    xsection += Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
  } 
  xsection *= millibarn; // parametrised in mb
  return xsection;
}

////////////////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4GlauberGribovCrossSection::GetNucleusRadius( const G4DynamicParticle* , 
                                               const G4Element* anElement)
{
  G4double At       = anElement->GetN();
  G4double oneThird = 1.0/3.0;
  G4double cubicrAt = std::pow (At, oneThird); 


  G4double R;  // = fRadiusConst*cubicrAt;
  /*  
  G4double tmp = std::pow( cubicrAt-1., 3.);
  tmp         += At;
  tmp         *= 0.5;

  if (At > 20.)   // 20.
  {
    R = fRadiusConst*std::pow (tmp, oneThird); 
  }
  else
  {
    R = fRadiusConst*cubicrAt; 
  }
  */
  
  R = fRadiusConst*cubicrAt;

  G4double meanA  = 21.;

  G4double tauA1  = 40.; 
  G4double tauA2  = 10.; 
  G4double tauA3  = 5.; 

  G4double a1 = 0.85;
  G4double b1 = 1. - a1;

  G4double b2 = 0.3;
  G4double b3 = 4.;

  if (At > 20.)   // 20.
  {
    R *= ( a1 + b1*std::exp( -(At - meanA)/tauA1) ); 
  }
  else if (At > 3.5)
  {
    R *= ( 1.0 + b2*( 1. - std::exp( (At - meanA)/tauA2) ) ); 
  }
  else 
  {
    R *= ( 1.0 + b3*( 1. - std::exp( (At - meanA)/tauA3) ) ); 
  }
  
  return R;
}
////////////////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4GlauberGribovCrossSection::GetNucleusRadius(G4double At)
{
  G4double oneThird = 1.0/3.0;
  G4double cubicrAt = std::pow (At, oneThird); 


  G4double R;  // = fRadiusConst*cubicrAt;

  /*
  G4double tmp = std::pow( cubicrAt-1., 3.);
  tmp         += At;
  tmp         *= 0.5;

  if (At > 20.)
  {
    R = fRadiusConst*std::pow (tmp, oneThird); 
  }
  else
  {
    R = fRadiusConst*cubicrAt; 
  }
  */

  R = fRadiusConst*cubicrAt;

  G4double meanA = 20.;
  G4double tauA  = 20.; 

  if (At > 20.)   // 20.
  {
    R *= ( 0.8 + 0.2*std::exp( -(At - meanA)/tauA) ); 
  }
  else
  {
    R *= ( 1.0 + 0.1*( 1. - std::exp( (At - meanA)/tauA) ) ); 
  }

  return R;
}

////////////////////////////////////////////////////////////////////////////////////
//
//

G4double G4GlauberGribovCrossSection::CalculateEcmValue( const G4double mp , 
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

G4double G4GlauberGribovCrossSection::CalcMandelstamS( const G4double mp , 
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
