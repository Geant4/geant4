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
  fLowerLimit( 1 * GeV ),
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

G4bool 
G4GlauberGribovCrossSection::IsApplicable(const G4DynamicParticle* aDP, const G4Element*)
{
  G4bool applicable  = false;
  G4int baryonNumber     = aDP->GetDefinition()->GetBaryonNumber();
  G4double kineticEnergy = aDP->GetKineticEnergy();

  const G4ParticleDefinition* theParticle = aDP->GetDefinition();
 
  if ( kineticEnergy / baryonNumber <= fUpperLimit &&
       ( theParticle == theProton    ||
         theParticle == theAProton   ||
         theParticle == theGamma     ||
         theParticle == thePiPlus    ||
         theParticle == thePiMinus   ||
         theParticle == theKPlus     ||
         theParticle == theKMinus    || 
         theParticle == theNeutron      ) ) applicable = true;
  return applicable;
}


////////////////////////////////////////////////////////////////////////////////////////
//
// Calculates total and elastic Xsc, derives inelatic as total - elastic accordong to
// Glauber model with Gribov correction calculated in the dipole approximation on
// light cone. Gaussian density helps to calculate rest integrals of the model.
// [1] B.Z. Kopeliovich, nucl-th/0306044 



G4double G4GlauberGribovCrossSection::
GetCrossSection(const G4DynamicParticle* aParticle, const G4Element* anElement, G4double )
{
  G4double xsection;
  G4double At           = anElement->GetN();
  G4double one_third = 1.0 / 3.0;
  G4double cubicrAt  = std::pow ( At , G4double(one_third) ); 

  G4double sigma     = GetHadronNucleaonXsc(aParticle, anElement);

  G4double R             = fRadiusConst*cubicrAt;
  G4double nucleusSquare = 2.*pi*R*R; 
  G4double ratio = sigma/nucleusSquare;

  xsection =  nucleusSquare*std::log( 1. + ratio );

  fTotalXsc = xsection;
   
  fElasticXsc = 0.5*( xsection - nucleusSquare*ratio/(1.+ratio) );

  if (fElasticXsc < 0.) fElasticXsc = 0.;

  fInelasticXsc = fTotalXsc - fElasticXsc;

  if (fInelasticXsc < 0.) fInelasticXsc = 0.;
  
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
  G4double xsection;

  G4double At = anElement->GetN();
  G4double Zt = anElement->GetZ();  


  G4double targ_mass = G4ParticleTable::GetParticleTable()->
  GetIonTable()->GetIonMass( G4int(Zt+0.5) , G4int(At+0.5) );
  G4double proj_mass     = aParticle->GetMass();
  G4double proj_momentum = aParticle->GetMomentum().mag();
  G4double sMand = CalcMandelstamS ( proj_mass , targ_mass , proj_momentum );

  sMand /= GeV*GeV;  // in GeV for parametrisation

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
