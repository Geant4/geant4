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
: fUpperLimit( 1000 * GeV ),
  fLowerLimit( 1 * GeV ),
  fRadiusConst( 1.1*fermi )
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

////////////////////////////////////////////////////////////////////////////////////////
//
//



G4double G4GlauberGribovCrossSection::
GetCrossSection(const G4DynamicParticle* aParticle, const G4Element* anElement, G4double )
{
  G4double xsection = 0.0;

  G4int Ap = aParticle->GetDefinition()->GetBaryonNumber();
  // G4int Zp = int ( aParticle->GetDefinition()->GetPDGCharge() / eplus + 0.5 ); 
  G4double ke_per_N = aParticle->GetKineticEnergy() / Ap; 

  G4int At = G4int ( anElement->GetN() + 0.5 );
  // G4int Zt = G4int ( anElement->GetZ() + 0.5 );  

// Apply energy check, if less than lower limit then 0 value is returned

  if (  ke_per_N < fLowerLimit )  return xsection;
 
  G4double one_third = 1.0 / 3.0;

  G4double cubicrAt = std::pow ( G4double(At) , G4double(one_third) ); 

  G4double sigma = GetHadronNucleaonXsc(aParticle, anElement);


  G4double R = fRadiusConst * cubicrAt;
  G4double nucleusSquare = pi*R*R; 

  xsection =  nucleusSquare* std::log( 1 + At*sigma/nucleusSquare );
   
  
  return xsection; 
}

/////////////////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4GlauberGribovCrossSection::GetHadronNucleaonXsc(const G4DynamicParticle* aParticle, 
                                                  const G4Element* anElement          )
{
  G4double xsection;

  G4int At = G4int ( anElement->GetN() + 0.5 );
  G4int Zt = G4int ( anElement->GetZ() + 0.5 );  


  G4double targ_mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass( Zt , At );
  G4double proj_mass     = aParticle->GetMass();
  G4double proj_momentum = aParticle->GetMomentum().mag();
  G4double Ecm = CalculateEcmValue ( proj_mass , targ_mass , proj_momentum );

  Ecm /= GeV;  // in GeV for parametrisation

  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  

  if(theParticle == theGamma) 
  {
    xsection = 0.0677*std::pow(Ecm,0.0808) + 0.129*std::pow(Ecm,-0.4525);
  } 
  else if(theParticle == theProton) 
  {
    xsection = 21.70*std::pow(Ecm,0.0808) + 56.08*std::pow(Ecm,-0.4525);
  } 
  else if(theParticle == theAProton) 
  {
    xsection = 21.70*std::pow(Ecm,0.0808) + 98.39*std::pow(Ecm,-0.4525);
  } 
  else if(theParticle == theNeutron) 
  {
    xsection = 21.70*std::pow(Ecm,0.0808) + 56.08*std::pow(Ecm,-0.4525); // ??? as proton ???
  } 
  else if(theParticle == thePiPlus) 
  {
    xsection = 13.63*std::pow(Ecm,0.0808) + 27.56*std::pow(Ecm,-0.4525);
  } 
  else if(theParticle == thePiMinus) 
  {
    xsection = 13.63*std::pow(Ecm,0.0808) + 36.02*std::pow(Ecm,-0.4525);
  } 
  else if(theParticle == theKPlus) 
  {
    xsection = 11.82*std::pow(Ecm,0.0808) + 8.15*std::pow(Ecm,-0.4525);
  } 
  else if(theParticle == theKMinus) 
  {
    xsection = 11.82*std::pow(Ecm,0.0808) + 26.36*std::pow(Ecm,-0.4525);
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
  G4double Pcm  = Plab * mt / Ecm;
  G4double KEcm = std::sqrt ( Pcm * Pcm + mp * mp ) - mp;

  return KEcm;
}

///////////////////////////////////////////////////////////////////////////////////////
//
//

G4double G4GlauberGribovCrossSection::CalculateCeValue( const G4double ke )
{
  // Calculate c value 
  // This value is indepenent from projectile and target particle 
  // ke is projectile kinetic energy per nucleon in the Lab system with MeV unit 
  // fitting function is made by T. Koi 
  // There are no data below 30 MeV/n in Kox et al., 

  G4double Ce; 
  G4double log10_ke = std::log10 ( ke );
   
  if ( log10_ke > 1.5 ) 
  {
      Ce = - 10.0 / std::pow ( G4double(log10_ke) , G4double(5) ) + 2.0;
  }
  else
  {
      Ce = ( - 10.0 / std::pow ( G4double(1.5) , G4double(5) ) + 2.0 ) / 
           std::pow ( G4double(1.5) , G4double(3) ) * 
           std::pow ( G4double(log10_ke) , G4double(3) );
  }
  return Ce;
}


//
//
///////////////////////////////////////////////////////////////////////////////////////
