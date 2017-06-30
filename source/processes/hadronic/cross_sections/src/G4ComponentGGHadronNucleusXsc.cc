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
// author: V. Grichine
// 
// 25.04.12 V. Grichine - first implementation

#include "G4ComponentGGHadronNucleusXsc.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4HadronNucleonXsc.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"

//////////////////////////////////////////////////////////////////////////////
//

G4ComponentGGHadronNucleusXsc::G4ComponentGGHadronNucleusXsc() 
 : G4VComponentCrossSection(Default_Name()),
//   fUpperLimit(100000*GeV),
   fLowerLimit(10.*MeV),// fLowerLimit(3*GeV),
   fRadiusConst(1.08*fermi),  // 1.1, 1.3 ?
   fTotalXsc(0.0), fElasticXsc(0.0), fInelasticXsc(0.0), fProductionXsc(0.0),
   fDiffractionXsc(0.0), fAxsc2piR2(0.0),fModelInLog(0.0)
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

  hnXsc = new G4HadronNucleonXsc();
}

///////////////////////////////////////////////////////////////////////////////////////
//
//

G4ComponentGGHadronNucleusXsc::~G4ComponentGGHadronNucleusXsc()
{
  if (hnXsc) delete hnXsc;
}

////////////////////////////////////////////////////////////////////

G4double G4ComponentGGHadronNucleusXsc::GetTotalIsotopeCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy,
				       G4int Z, G4int A)
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fTotalXsc = GetIsoCrossSection(aDP, Z, A);
  delete aDP;

  return fTotalXsc;
}

//////////////////////////////////////////////////////////////////////

G4double G4ComponentGGHadronNucleusXsc::GetTotalElementCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy, 
				       G4int Z, G4double A)
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fTotalXsc = GetIsoCrossSection(aDP, Z, G4int(A));
  delete aDP;

  return fTotalXsc;
}

////////////////////////////////////////////////////////////////////

G4double G4ComponentGGHadronNucleusXsc::GetInelasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4int A)
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fTotalXsc = GetIsoCrossSection(aDP, Z, A);
  delete aDP;

  return fInelasticXsc;
}

////////////////////////////////////////////////////////////////////

G4double G4ComponentGGHadronNucleusXsc::GetProductionIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4int A)
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fTotalXsc = GetIsoCrossSection(aDP, Z, A);
  delete aDP;

  return fProductionXsc;
}

/////////////////////////////////////////////////////////////////////

G4double G4ComponentGGHadronNucleusXsc::GetInelasticElementCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4double A)
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fTotalXsc = GetIsoCrossSection(aDP, Z, G4int(A));
  delete aDP;

  return fInelasticXsc;
}

/////////////////////////////////////////////////////////////////////

G4double G4ComponentGGHadronNucleusXsc::GetProductionElementCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4double A)
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fTotalXsc = GetIsoCrossSection(aDP, Z, G4int(A));
  delete aDP;

  return fProductionXsc;
}

//////////////////////////////////////////////////////////////////

G4double G4ComponentGGHadronNucleusXsc::GetElasticElementCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4double A)
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fTotalXsc = GetIsoCrossSection(aDP, Z, G4int(A));
  delete aDP;

  return fElasticXsc;
}

///////////////////////////////////////////////////////////////////

G4double G4ComponentGGHadronNucleusXsc::GetElasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4int A)
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fTotalXsc = GetIsoCrossSection(aDP, Z, A);
  delete aDP;

  return fElasticXsc;
}

////////////////////////////////////////////////////////////////
 
G4double G4ComponentGGHadronNucleusXsc::ComputeQuasiElasticRatio(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4int A)
{
  G4DynamicParticle* aDP = new G4DynamicParticle(aParticle,G4ParticleMomentum(1.,0.,0.), 
                                                kinEnergy);
  fTotalXsc = GetIsoCrossSection(aDP, Z, A);
  delete aDP;
  G4double ratio = 0.;

  if(fInelasticXsc > 0.)
  {
    ratio = (fInelasticXsc - fProductionXsc)/fInelasticXsc;
    if(ratio < 0.) ratio = 0.;
  }
  return ratio;
}
 



////////////////////////////////////////////////////////////////////////////////////////

G4bool 
G4ComponentGGHadronNucleusXsc::IsIsoApplicable(const G4DynamicParticle* aDP, 
					     G4int Z, G4int /*A*/, 
					     const G4Element*,
					     const G4Material*)
{
  G4bool applicable      = false;
  // G4int baryonNumber     = aDP->GetDefinition()->GetBaryonNumber();
  G4double kineticEnergy = aDP->GetKineticEnergy();

  const G4ParticleDefinition* theParticle = aDP->GetDefinition();
 
  if ( 
       Z >= 1  // >=  H for kaons
       &&     
       ( 
         kineticEnergy  >= fLowerLimit 
         &&
	 //         Z > 1 &&      // >=  He
         ( 
           theParticle == theAProton   ||
           theParticle == theGamma     ||
           theParticle == theSMinus    ||  
           theParticle == theProton    ||
           theParticle == theNeutron   ||   
           theParticle == thePiPlus    ||
           theParticle == thePiMinus       
         )
       )  
     ) 
     applicable = true;

  if ( 
       Z >= 1  // >=  H for kaons
       &&
       ( 
         kineticEnergy  >= 0.01*fLowerLimit 
         &&
         ( 
           theParticle == theKPlus     ||
           theParticle == theKMinus    || 
           theParticle == theK0L       ||
           theParticle == theK0S       
         )    
       )    
     ) 
     applicable = true;

  return applicable;
}

////////////////////////////////////////////////////////////////////////////////////////
//
// Calculates total and inelastic Xsc, derives elastic as total - inelastic accordong to
// Glauber model with Gribov correction calculated in the dipole approximation on
// light cone. Gaussian density of point-like nucleons helps to calculate rest integrals of the model.
// [1] B.Z. Kopeliovich, nucl-th/0306044 + simplification above

G4double 
G4ComponentGGHadronNucleusXsc::GetIsoCrossSection(const G4DynamicParticle* aParticle, 
						G4int Z, G4int A,  
						const G4Isotope*,
						const G4Element*,
						const G4Material*)
{
  G4double xsection, sigma, cofInelastic, cofTotal, nucleusSquare, ratio;
  G4double hpInXsc(0.), hnInXsc(0.);
  G4double R             = GetNucleusRadius(A); 
  
  G4int N = A - Z;              // number of neutrons
  if (N < 0) N = 0;

  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();

  if( theParticle == theProton  || 
      theParticle == theNeutron ||
      theParticle == thePiPlus  || 
      theParticle == thePiMinus      )
  {
    // sigma        = GetHadronNucleonXscNS(aParticle, A, Z);

    sigma = Z*hnXsc->GetHadronNucleonXscNS(aParticle, theProton);

    hpInXsc = hnXsc->GetInelasticHadronNucleonXsc();

    sigma += N*hnXsc->GetHadronNucleonXscNS(aParticle, theNeutron);

    hnInXsc = hnXsc->GetInelasticHadronNucleonXsc();

    cofInelastic = 2.4;
    cofTotal     = 2.0;
  }
  else if( theParticle == theKPlus   || 
           theParticle == theKMinus  || 
           theParticle == theK0S     || 
           theParticle == theK0L        ) 
  {
    // sigma        = GetKaonNucleonXscVector(aParticle, A, Z);

    sigma = Z*hnXsc->GetKaonNucleonXscGG(aParticle, theProton);

    hpInXsc = hnXsc->GetInelasticHadronNucleonXsc();

    sigma += N*hnXsc->GetKaonNucleonXscGG(aParticle, theNeutron);

    hnInXsc = hnXsc->GetInelasticHadronNucleonXsc();

    cofInelastic = 2.2;
    cofTotal     = 2.0;
    R = 1.3*fermi;
    R *= G4Pow::GetInstance()->powA(G4double(A), 0.3333);
  }
  else
  {
    sigma        = GetHadronNucleonXscNS(aParticle, A, Z);
    cofInelastic = 2.2;
    cofTotal     = 2.0;
  }
  // cofInelastic = 2.0;

  if( A > 1 )
  { 
    nucleusSquare = cofTotal*pi*R*R;   // basically 2piRR
    ratio = sigma/nucleusSquare;

    xsection =  nucleusSquare*G4Log( 1. + ratio );

    xsection *= GetParticleBarCorTot(theParticle, Z);

    fTotalXsc = xsection;

    // inelastic xsc

    fAxsc2piR2 = cofInelastic*ratio;

    fModelInLog = G4Log( 1. + fAxsc2piR2 );

    fInelasticXsc = nucleusSquare*fModelInLog/cofInelastic;

    fInelasticXsc *= GetParticleBarCorIn(theParticle, Z);

    fElasticXsc   = fTotalXsc - fInelasticXsc;

    if( fElasticXsc < 0. ) fElasticXsc = 0.;
    
    G4double difratio = ratio/(1.+ratio);

    fDiffractionXsc = 0.5*nucleusSquare*( difratio - G4Log( 1. + difratio ) );


    // sigma = GetHNinelasticXsc(aParticle, A, Z);

    sigma = Z*hpInXsc + N*hnInXsc;

    ratio = sigma/nucleusSquare;

    fProductionXsc = nucleusSquare*G4Log( 1. + cofInelastic*ratio )/cofInelastic;

    fProductionXsc *= GetParticleBarCorIn(theParticle, Z);

    if (fElasticXsc < 0.) fElasticXsc = 0.;
  }
  else // H
  {
    fTotalXsc = sigma;
    xsection  = sigma;

    fInelasticXsc = hnXsc->GetInelasticHadronNucleonXsc();

    if ( theParticle != theAProton ) 
    {
     fElasticXsc = hnXsc->GetElasticHadronNucleonXsc();

     //      sigma         = GetHNinelasticXsc(aParticle, A, Z);
     // fInelasticXsc = sigma;
     // fElasticXsc   = fTotalXsc - fInelasticXsc;      
    }
    else if( theParticle == theKPlus || 
             theParticle == theKMinus  || 
             theParticle == theK0S     || 
             theParticle == theK0L        ) 
    { 
      fInelasticXsc = hpInXsc;
      fElasticXsc   = fTotalXsc - fInelasticXsc;
    }   
    else
    {
      fInelasticXsc = hpInXsc;
      fElasticXsc   = fTotalXsc - fInelasticXsc;
    }
    if (fElasticXsc < 0.) fElasticXsc = 0.;
      
  }
  return xsection; 
}

//////////////////////////////////////////////////////////////////////////
//
// Return single-diffraction/inelastic cross-section ratio

G4double G4ComponentGGHadronNucleusXsc::
GetRatioSD(const G4DynamicParticle* aParticle, G4int A, G4int Z)
{
  G4double sigma, cofInelastic, cofTotal, nucleusSquare, ratio;
  G4double R             = GetNucleusRadius(A); 

  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();

  if( theParticle == theProton  || 
      theParticle == theNeutron ||
      theParticle == thePiPlus  || 
      theParticle == thePiMinus      )
  {
    sigma        = GetHadronNucleonXscNS(aParticle, A, Z);
    cofInelastic = 2.4;
    cofTotal     = 2.0;
  }
  else
  {
    sigma        = GetHadronNucleonXscNS(aParticle, A, Z);
    cofInelastic = 2.2;
    cofTotal     = 2.0;
  }
  nucleusSquare = cofTotal*pi*R*R;   // basically 2piRR
  ratio = sigma/nucleusSquare;

  fInelasticXsc = nucleusSquare*G4Log( 1. + cofInelastic*ratio )/cofInelastic;
   
  G4double difratio = ratio/(1.+ratio);

  fDiffractionXsc = 0.5*nucleusSquare*( difratio - G4Log( 1. + difratio ) );

  if (fInelasticXsc > 0.) ratio = fDiffractionXsc/fInelasticXsc;
  else                    ratio = 0.;

  return ratio; 
}

//////////////////////////////////////////////////////////////////////////
//
// Return suasi-elastic/inelastic cross-section ratio

G4double G4ComponentGGHadronNucleusXsc::
GetRatioQE(const G4DynamicParticle* aParticle, G4int A, G4int Z)
{
  G4double sigma, cofInelastic, cofTotal, nucleusSquare, ratio;
  G4double R             = GetNucleusRadius(A); 

  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();

  if( theParticle == theProton  || 
      theParticle == theNeutron ||
      theParticle == thePiPlus  || 
      theParticle == thePiMinus      )
  {
    sigma        = GetHadronNucleonXscNS(aParticle, A, Z);
    cofInelastic = 2.4;
    cofTotal     = 2.0;
  }
  else
  {
    sigma        = GetHadronNucleonXscNS(aParticle, A, Z);
    cofInelastic = 2.2;
    cofTotal     = 2.0;
  }
  nucleusSquare = cofTotal*pi*R*R;   // basically 2piRR
  ratio = sigma/nucleusSquare;

  fInelasticXsc = nucleusSquare*G4Log( 1. + cofInelastic*ratio )/cofInelastic;

  sigma = GetHNinelasticXsc(aParticle, A, Z);
  ratio = sigma/nucleusSquare;

  fProductionXsc = nucleusSquare*G4Log( 1. + cofInelastic*ratio )/cofInelastic;

  if (fInelasticXsc > fProductionXsc) ratio = (fInelasticXsc-fProductionXsc)/fInelasticXsc;
  else                                ratio = 0.;
  if ( ratio < 0. )                   ratio = 0.;

  return ratio; 
}

/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to differnt parametrisations:
// [2] E. Levin, hep-ph/9710546
// [3] U. Dersch, et al, hep-ex/9910052
// [4] M.J. Longo, et al, Phys.Rev.Lett. 33 (1974) 725 

G4double 
G4ComponentGGHadronNucleusXsc::GetHadronNucleonXsc(const G4DynamicParticle* aParticle, 
                                                 const G4Element* anElement)
{
  G4int At = G4lrint(anElement->GetN());  // number of nucleons 
  G4int Zt = G4lrint(anElement->GetZ());  // number of protons

  return GetHadronNucleonXsc(aParticle, At, Zt);
}

/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to differnt parametrisations:
// [2] E. Levin, hep-ph/9710546
// [3] U. Dersch, et al, hep-ex/9910052
// [4] M.J. Longo, et al, Phys.Rev.Lett. 33 (1974) 725 

G4double 
G4ComponentGGHadronNucleusXsc::GetHadronNucleonXsc(const G4DynamicParticle* aParticle, 
                                                 G4int At, G4int /*Zt*/)
{
  G4double xsection;

  //G4double targ_mass = G4NucleiProperties::GetNuclearMass(At, Zt);

  G4double targ_mass = 0.939*GeV;  // ~mean neutron and proton ???

  G4double proj_mass     = aParticle->GetMass();
  G4double proj_momentum = aParticle->GetMomentum().mag();
  G4double sMand = CalcMandelstamS ( proj_mass , targ_mass , proj_momentum );

  sMand /= GeV*GeV;  // in GeV for parametrisation
  proj_momentum /= GeV;

  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  
  G4double aa = At;

  if(theParticle == theGamma) 
  {
    xsection = aa*(0.0677*G4Pow::GetInstance()->powA(sMand,0.0808) + 0.129*G4Pow::GetInstance()->powA(sMand,-0.4525));
  } 
  else if(theParticle == theNeutron) // as proton ??? 
  {
    xsection = aa*(21.70*G4Pow::GetInstance()->powA(sMand,0.0808) + 56.08*G4Pow::GetInstance()->powA(sMand,-0.4525));
  } 
  else if(theParticle == theProton) 
  {
    xsection = aa*(21.70*G4Pow::GetInstance()->powA(sMand,0.0808) + 56.08*G4Pow::GetInstance()->powA(sMand,-0.4525));
    // xsection = At*( 49.51*G4Pow::GetInstance()->powA(sMand,-0.097) + 0.314*G4Log(sMand)*G4Log(sMand) );
    // xsection = At*( 38.4 + 0.85*std::abs(G4Pow::GetInstance()->powA(log(sMand),1.47)) );
  } 
  else if(theParticle == theAProton) 
  {
    xsection = aa*( 21.70*G4Pow::GetInstance()->powA(sMand,0.0808) + 98.39*G4Pow::GetInstance()->powA(sMand,-0.4525));
  } 
  else if(theParticle == thePiPlus) 
  {
    xsection = aa*(13.63*G4Pow::GetInstance()->powA(sMand,0.0808) + 27.56*G4Pow::GetInstance()->powA(sMand,-0.4525));
  } 
  else if(theParticle == thePiMinus) 
  {
    // xsection = At*( 55.2*G4Pow::GetInstance()->powA(sMand,-0.255) + 0.346*G4Log(sMand)*G4Log(sMand) );
    xsection = aa*(13.63*G4Pow::GetInstance()->powA(sMand,0.0808) + 36.02*G4Pow::GetInstance()->powA(sMand,-0.4525));
  } 
  else if(theParticle == theKPlus) 
  {
    xsection = aa*(11.82*G4Pow::GetInstance()->powA(sMand,0.0808) + 8.15*G4Pow::GetInstance()->powA(sMand,-0.4525));
  } 
  else if(theParticle == theKMinus) 
  {
    xsection = aa*(11.82*G4Pow::GetInstance()->powA(sMand,0.0808) + 26.36*G4Pow::GetInstance()->powA(sMand,-0.4525));
  }
  else  // as proton ??? 
  {
    xsection = aa*(21.70*G4Pow::GetInstance()->powA(sMand,0.0808) + 56.08*G4Pow::GetInstance()->powA(sMand,-0.4525));
  } 
  xsection *= millibarn;
  return xsection;
}


/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to PDG parametrisation (2005):
// http://pdg.lbl.gov/2006/reviews/hadronicrpp.pdf

G4double 
G4ComponentGGHadronNucleusXsc::GetHadronNucleonXscPDG(const G4DynamicParticle* aParticle, 
                                                    const G4Element* anElement)
{
  G4int At = G4lrint(anElement->GetN());  // number of nucleons 
  G4int Zt = G4lrint(anElement->GetZ());  // number of protons

  return GetHadronNucleonXscPDG(aParticle, At, Zt);
}




/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to PDG parametrisation (2005):
// http://pdg.lbl.gov/2006/reviews/hadronicrpp.pdf
//  At = number of nucleons,  Zt = number of protons 

G4double 
G4ComponentGGHadronNucleusXsc::GetHadronNucleonXscPDG(const G4DynamicParticle* aParticle, 
                                                    G4int At, G4int Zt)
{
  G4double xsection;

  G4int Nt = At-Zt;              // number of neutrons
  if (Nt < 0) Nt = 0;
  
  G4double zz = Zt;
  G4double aa = At;
  G4double nn = Nt;

  G4double targ_mass = G4ParticleTable::GetParticleTable()->
    GetIonTable()->GetIonMass(Zt, At);

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
    xsection = zz*( 35.80 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 40.15*G4Pow::GetInstance()->powA(sMand,-eta1) - 30.*G4Pow::GetInstance()->powA(sMand,-eta2));
    xsection  += nn*( 35.45 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
		      + 42.53*G4Pow::GetInstance()->powA(sMand,-eta1) - 33.34*G4Pow::GetInstance()->powA(sMand,-eta2)); // pp for nn
  } 
  else if(theParticle == theProton) 
  {
      
      xsection  = zz*( 35.45 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 42.53*G4Pow::GetInstance()->powA(sMand,-eta1) - 33.34*G4Pow::GetInstance()->powA(sMand,-eta2));

      xsection += nn*( 35.80 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 40.15*G4Pow::GetInstance()->powA(sMand,-eta1) - 30.*G4Pow::GetInstance()->powA(sMand,-eta2));
  } 
  else if(theParticle == theAProton) 
  {
    xsection  = zz*( 35.45 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 42.53*G4Pow::GetInstance()->powA(sMand,-eta1) + 33.34*G4Pow::GetInstance()->powA(sMand,-eta2));

    xsection += nn*( 35.80 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 40.15*G4Pow::GetInstance()->powA(sMand,-eta1) + 30.*G4Pow::GetInstance()->powA(sMand,-eta2));
  } 
  else if(theParticle == thePiPlus) 
  {
    xsection  = aa*( 20.86 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 19.24*G4Pow::GetInstance()->powA(sMand,-eta1) - 6.03*G4Pow::GetInstance()->powA(sMand,-eta2));
  } 
  else if(theParticle == thePiMinus) 
  {
    xsection  = aa*( 20.86 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 19.24*G4Pow::GetInstance()->powA(sMand,-eta1) + 6.03*G4Pow::GetInstance()->powA(sMand,-eta2));
  } 
  else if(theParticle == theKPlus || theParticle == theK0L ) 
  {
    xsection  = zz*( 17.91 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 7.14*G4Pow::GetInstance()->powA(sMand,-eta1) - 13.45*G4Pow::GetInstance()->powA(sMand,-eta2));

    xsection += nn*( 17.87 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 5.17*G4Pow::GetInstance()->powA(sMand,-eta1) - 7.23*G4Pow::GetInstance()->powA(sMand,-eta2));
  } 
  else if(theParticle == theKMinus || theParticle == theK0S ) 
  {
    xsection  = zz*( 17.91 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 7.14*G4Pow::GetInstance()->powA(sMand,-eta1) + 13.45*G4Pow::GetInstance()->powA(sMand,-eta2));

    xsection += nn*( 17.87 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 5.17*G4Pow::GetInstance()->powA(sMand,-eta1) + 7.23*G4Pow::GetInstance()->powA(sMand,-eta2));
  }
  else if(theParticle == theSMinus) 
  {
    xsection  = aa*( 35.20 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          - 199.*G4Pow::GetInstance()->powA(sMand,-eta1) + 264.*G4Pow::GetInstance()->powA(sMand,-eta2));
  } 
  else if(theParticle == theGamma) // modify later on
  {
    xsection  = aa*( 0.0 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 0.032*G4Pow::GetInstance()->powA(sMand,-eta1) - 0.0*G4Pow::GetInstance()->powA(sMand,-eta2));
   
  } 
  else  // as proton ??? 
  {
    xsection  = zz*( 35.45 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 42.53*G4Pow::GetInstance()->powA(sMand,-eta1) - 33.34*G4Pow::GetInstance()->powA(sMand,-eta2));

    xsection += nn*( 35.80 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 40.15*G4Pow::GetInstance()->powA(sMand,-eta1) - 30.*G4Pow::GetInstance()->powA(sMand,-eta2));
  } 
  xsection *= millibarn; // parametrised in mb
  return xsection;
}


/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon cross-section based on N. Starkov parametrisation of
// data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database

G4double 
G4ComponentGGHadronNucleusXsc::GetHadronNucleonXscNS(const G4DynamicParticle* aParticle, 
                                                   const G4Element* anElement)
{
  G4int At = G4lrint(anElement->GetN());  // number of nucleons 
  G4int Zt = G4lrint(anElement->GetZ());  // number of protons

  return GetHadronNucleonXscNS(aParticle, At, Zt);
}




/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon cross-section based on N. Starkov parametrisation of
// data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database

G4double 
G4ComponentGGHadronNucleusXsc::GetHadronNucleonXscNS(const G4DynamicParticle* aParticle, 
                                                   G4int At, G4int Zt)
{
  G4double xsection(0);
  // G4double Delta;   DHW 19 May 2011: variable set but not used
  G4double A0, B0;
  G4double hpXscv(0);
  G4double hnXscv(0);

  G4int Nt = At-Zt;              // number of neutrons
  if (Nt < 0) Nt = 0;  

  G4double aa = At;
  G4double zz = Zt;
  G4double nn = Nt;

  G4double targ_mass = G4ParticleTable::GetParticleTable()->
  GetIonTable()->GetIonMass(Zt, At);

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
    if( proj_momentum >= 373.)
    {
      return GetHadronNucleonXscPDG(aParticle,At,Zt);
    }
    else if( proj_momentum >= 10.)
    // if( proj_momentum >= 2.)
    {
      //  Delta = 1.;  // DHW 19 May 2011: variable set but not used
      // if( proj_energy < 40. ) Delta = 0.916+0.0021*proj_energy;

      //AR-12Aug2016  if(proj_momentum >= 10.)
      {
        B0 = 7.5;
        A0 = 100. - B0*G4Log(3.0e7);

        xsection = A0 + B0*G4Log(proj_energy) - 11
                  + 103*G4Pow::GetInstance()->powA(2*0.93827*proj_energy + proj_mass*proj_mass+
                     0.93827*0.93827,-0.165);        //  mb
      }
      xsection *= zz + nn;
    }
    else
    {
      // nn to be pp

      if( proj_momentum < 0.73 )
      {
        hnXscv = 23 + 50*( G4Pow::GetInstance()->powA( G4Log(0.73/proj_momentum), 3.5 ) );
      }
      else if( proj_momentum < 1.05  )
      {
       hnXscv = 23 + 40*(G4Log(proj_momentum/0.73))*
                         (G4Log(proj_momentum/0.73));
      }
      else  // if( proj_momentum < 10.  )
      {
         hnXscv = 39.0+
              75*(proj_momentum - 1.2)/(G4Pow::GetInstance()->powA(proj_momentum,3.0) + 0.15);
      }
      // pn to be np

      if( proj_momentum < 0.8 )
      {
        hpXscv = 33+30*G4Pow::GetInstance()->powA(G4Log(proj_momentum/1.3),4.0);
      }      
      else if( proj_momentum < 1.4 )
      {
        hpXscv = 33+30*G4Pow::GetInstance()->powA(G4Log(proj_momentum/0.95),2.0);
      }
      else    // if( proj_momentum < 10.  )
      {
        hpXscv = 33.3+
              20.8*(G4Pow::GetInstance()->powA(proj_momentum,2.0)-1.35)/
                 (G4Pow::GetInstance()->powA(proj_momentum,2.50)+0.95);
      }
      xsection = hpXscv*zz + hnXscv*nn;
    }
  } 
  else if(theParticle == theProton) 
  {
    if( proj_momentum >= 373.)
    {
      return GetHadronNucleonXscPDG(aParticle,At,Zt);
    }
    else if( proj_momentum >= 10.)
    // if( proj_momentum >= 2.)
    {
      // Delta = 1.;  DHW 19 May 2011: variable set but not used
      // if( proj_energy < 40. ) Delta = 0.916+0.0021*proj_energy;

      //AR-12Aug2016  if(proj_momentum >= 10.)
      {
        B0 = 7.5;
        A0 = 100. - B0*G4Log(3.0e7);

        xsection = A0 + B0*G4Log(proj_energy) - 11
                  + 103*G4Pow::GetInstance()->powA(2*0.93827*proj_energy + proj_mass*proj_mass+
                     0.93827*0.93827,-0.165);        //  mb
      }
      xsection *= zz + nn;
    }
    else
    {
      // pp

      if( proj_momentum < 0.73 )
      {
        hpXscv = 23 + 50*( G4Pow::GetInstance()->powA( G4Log(0.73/proj_momentum), 3.5 ) );
      }
      else if( proj_momentum < 1.05  )
      {
       hpXscv = 23 + 40*(G4Log(proj_momentum/0.73))*
                         (G4Log(proj_momentum/0.73));
      }
      else    // if( proj_momentum < 10.  )
      {
         hpXscv = 39.0+
              75*(proj_momentum - 1.2)/(G4Pow::GetInstance()->powA(proj_momentum,3.0) + 0.15);
      }
      // pn to be np

      if( proj_momentum < 0.8 )
      {
        hnXscv = 33+30*G4Pow::GetInstance()->powA(G4Log(proj_momentum/1.3),4.0);
      }      
      else if( proj_momentum < 1.4 )
      {
        hnXscv = 33+30*G4Pow::GetInstance()->powA(G4Log(proj_momentum/0.95),2.0);
      }
      else   // if( proj_momentum < 10.  )
      {
        hnXscv = 33.3+
              20.8*(G4Pow::GetInstance()->powA(proj_momentum,2.0)-1.35)/
                 (G4Pow::GetInstance()->powA(proj_momentum,2.50)+0.95);
      }
      xsection = hpXscv*zz + hnXscv*nn;
      // xsection = hpXscv*(Zt + Nt);
      // xsection = hnXscv*(Zt + Nt);
    }    
    // xsection *= 0.95;
  } 
  else if( theParticle == theAProton ) 
  {
    // xsection  = Zt*( 35.45 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
    //                       + 42.53*G4Pow::GetInstance()->powA(sMand,-eta1) + 33.34*G4Pow::GetInstance()->powA(sMand,-eta2));

    // xsection += Nt*( 35.80 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
    //                    + 40.15*G4Pow::GetInstance()->powA(sMand,-eta1) + 30.*G4Pow::GetInstance()->powA(sMand,-eta2));

    G4double logP = G4Log(proj_momentum);

    if( proj_momentum <= 1.0 )
    {
      xsection  = zz*(65.55 + 53.84/(proj_momentum+1.e-6)  );
    }
    else
    {
      xsection  = zz*( 41.1 + 77.2*G4Pow::GetInstance()->powA( proj_momentum, -0.68) 
                       + 0.293*logP*logP - 1.82*logP );
    }
    if ( nn > 0.)  
    {
      xsection += nn*( 41.9 + 96.2*G4Pow::GetInstance()->powA( proj_momentum, -0.99) - 0.154*logP);
    }
    else // H
    {
      fInelasticXsc =   38.0 + 38.0*G4Pow::GetInstance()->powA( proj_momentum, -0.96) 
	                - 0.169*logP*logP;
      fInelasticXsc *=  millibarn;
    }    
  } 
  else if( theParticle == thePiPlus ) 
  {
    if(proj_momentum < 0.4)
    {
      G4double Ex3 = 180*G4Exp(-(proj_momentum-0.29)*(proj_momentum-0.29)/0.085/0.085);
      hpXscv      = Ex3+20.0;
    }
    else if( proj_momentum < 1.15 )
    {
      G4double Ex4 = 88*(G4Log(proj_momentum/0.75))*(G4Log(proj_momentum/0.75));
      hpXscv = Ex4+14.0;
    }
    else if(proj_momentum < 3.5)
    {
      G4double Ex1 = 3.2*G4Exp(-(proj_momentum-2.55)*(proj_momentum-2.55)/0.55/0.55);
      G4double Ex2 = 12*G4Exp(-(proj_momentum-1.47)*(proj_momentum-1.47)/0.225/0.225);
      hpXscv = Ex1+Ex2+27.5;
    }
    else //  if(proj_momentum > 3.5) // mb
    {
      hpXscv = 10.6+2.*G4Log(proj_energy)+25*G4Pow::GetInstance()->powA(proj_energy,-0.43);
    }
    // pi+n = pi-p??

    if(proj_momentum < 0.37)
    {
      hnXscv = 28.0 + 40*G4Exp(-(proj_momentum-0.29)*(proj_momentum-0.29)/0.07/0.07);
    }
    else if(proj_momentum<0.65)
    {
       hnXscv = 26+110*(G4Log(proj_momentum/0.48))*(G4Log(proj_momentum/0.48));
    }
    else if(proj_momentum<1.3)
    {
      hnXscv = 36.1+
                10*G4Exp(-(proj_momentum-0.72)*(proj_momentum-0.72)/0.06/0.06)+
                24*G4Exp(-(proj_momentum-1.015)*(proj_momentum-1.015)/0.075/0.075);
    }
    else if(proj_momentum<3.0)
    {
      hnXscv = 36.1+0.079-4.313*G4Log(proj_momentum)+
                3*G4Exp(-(proj_momentum-2.1)*(proj_momentum-2.1)/0.4/0.4)+
                1.5*G4Exp(-(proj_momentum-1.4)*(proj_momentum-1.4)/0.12/0.12);
    }
    else   // mb
    {
      hnXscv = 10.6+2*G4Log(proj_energy)+30*G4Pow::GetInstance()->powA(proj_energy,-0.43); 
    }
    xsection = hpXscv*zz + hnXscv*nn;
  } 
  else if(theParticle == thePiMinus) 
  {
    // pi-n = pi+p??

    if(proj_momentum < 0.4)
    {
      G4double Ex3 = 180*G4Exp(-(proj_momentum-0.29)*(proj_momentum-0.29)/0.085/0.085);
      hnXscv      = Ex3+20.0;
    }
    else if(proj_momentum < 1.15)
    {
      G4double Ex4 = 88*(G4Log(proj_momentum/0.75))*(G4Log(proj_momentum/0.75));
      hnXscv = Ex4+14.0;
    }
    else if(proj_momentum < 3.5)
    {
      G4double Ex1 = 3.2*G4Exp(-(proj_momentum-2.55)*(proj_momentum-2.55)/0.55/0.55);
      G4double Ex2 = 12*G4Exp(-(proj_momentum-1.47)*(proj_momentum-1.47)/0.225/0.225);
      hnXscv = Ex1+Ex2+27.5;
    }
    else //  if(proj_momentum > 3.5) // mb
    {
      hnXscv = 10.6+2.*G4Log(proj_energy)+25*G4Pow::GetInstance()->powA(proj_energy,-0.43);
    }
    // pi-p

    if(proj_momentum < 0.37)
    {
      hpXscv = 28.0 + 40*G4Exp(-(proj_momentum-0.29)*(proj_momentum-0.29)/0.07/0.07);
    }
    else if(proj_momentum<0.65)
    {
       hpXscv = 26+110*(G4Log(proj_momentum/0.48))*(G4Log(proj_momentum/0.48));
    }
    else if(proj_momentum<1.3)
    {
      hpXscv = 36.1+
                10*G4Exp(-(proj_momentum-0.72)*(proj_momentum-0.72)/0.06/0.06)+
                24*G4Exp(-(proj_momentum-1.015)*(proj_momentum-1.015)/0.075/0.075);
    }
    else if(proj_momentum<3.0)
    {
      hpXscv = 36.1+0.079-4.313*G4Log(proj_momentum)+
                3*G4Exp(-(proj_momentum-2.1)*(proj_momentum-2.1)/0.4/0.4)+
                1.5*G4Exp(-(proj_momentum-1.4)*(proj_momentum-1.4)/0.12/0.12);
    }
    else   // mb
    {
      hpXscv = 10.6+2*G4Log(proj_energy)+30*G4Pow::GetInstance()->powA(proj_energy,-0.43); 
    }
    xsection = hpXscv*zz + hnXscv*nn;
  } 
  else if(theParticle == theKPlus) 
  {
    xsection  = zz*( 17.91 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 7.14*G4Pow::GetInstance()->powA(sMand,-eta1) - 13.45*G4Pow::GetInstance()->powA(sMand,-eta2));

    xsection += nn*( 17.87 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 5.17*G4Pow::GetInstance()->powA(sMand,-eta1) - 7.23*G4Pow::GetInstance()->powA(sMand,-eta2));
  } 
  else if(theParticle == theKMinus) 
  {
    xsection  = zz*( 17.91 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 7.14*G4Pow::GetInstance()->powA(sMand,-eta1) + 13.45*G4Pow::GetInstance()->powA(sMand,-eta2));

    xsection += nn*( 17.87 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 5.17*G4Pow::GetInstance()->powA(sMand,-eta1) + 7.23*G4Pow::GetInstance()->powA(sMand,-eta2));
  }
  else if(theParticle == theSMinus) 
  {
    xsection  = aa*( 35.20 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          - 199.*G4Pow::GetInstance()->powA(sMand,-eta1) + 264.*G4Pow::GetInstance()->powA(sMand,-eta2));
  } 
  else if(theParticle == theGamma) // modify later on
  {
    xsection  = aa*( 0.0 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 0.032*G4Pow::GetInstance()->powA(sMand,-eta1) - 0.0*G4Pow::GetInstance()->powA(sMand,-eta2));
   
  } 
  else  // as proton ??? 
  {
    xsection  = zz*( 35.45 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 42.53*G4Pow::GetInstance()->powA(sMand,-eta1) - 33.34*G4Pow::GetInstance()->powA(sMand,-eta2));

    xsection += nn*( 35.80 + B*G4Pow::GetInstance()->powA(G4Log(sMand/s0),2.) 
                          + 40.15*G4Pow::GetInstance()->powA(sMand,-eta1) - 30.*G4Pow::GetInstance()->powA(sMand,-eta2));
  } 
  xsection *= millibarn; // parametrised in mb
  return xsection;
}

/*
G4double 
G4ComponentGGHadronNucleusXsc::GetKaonNucleonXscVector(const G4DynamicParticle* aParticle, 
                                                   G4int At, G4int Zt)
{
  G4double Tkin, logTkin, xsc, xscP, xscN;
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();

  G4int Nt = At-Zt;              // number of neutrons
  if (Nt < 0) Nt = 0;  

  Tkin = aParticle->GetKineticEnergy(); // Tkin in MeV

  if( Tkin > 70*GeV ) return GetHadronNucleonXscPDG(aParticle,At,Zt);

  logTkin = G4Log(Tkin); // Tkin in MeV!!!

 if( theParticle == theKPlus )
 {
   xscP = hnXsc->GetKpProtonTotXscVector(logTkin);
   xscN = hnXsc->GetKpNeutronTotXscVector(logTkin);
 }
 else if( theParticle == theKMinus )
 {
   xscP = hnXsc->GetKmProtonTotXscVector(logTkin);
   xscN = hnXsc->GetKmNeutronTotXscVector(logTkin);
 }
 else // K-zero as half of K+ and K-
 {
   xscP = (hnXsc->GetKpProtonTotXscVector(logTkin)+hnXsc->GetKmProtonTotXscVector(logTkin))*0.5;
   xscN = (hnXsc->GetKpNeutronTotXscVector(logTkin)+hnXsc->GetKmNeutronTotXscVector(logTkin))*0.5;
 }
 xsc = xscP*Zt + xscN*Nt;
  return xsc;
}
*/

/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon inelastic cross-section based on proper parametrisation 

G4double 
G4ComponentGGHadronNucleusXsc::GetHNinelasticXsc(const G4DynamicParticle* aParticle, 
                                               const G4Element* anElement)
{
  G4int At = G4lrint(anElement->GetN());  // number of nucleons 
  G4int Zt = G4lrint(anElement->GetZ());  // number of protons

  return GetHNinelasticXsc(aParticle, At, Zt);
}

/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon inelastic cross-section based on FTF-parametrisation 

G4double 
G4ComponentGGHadronNucleusXsc::GetHNinelasticXsc(const G4DynamicParticle* aParticle, 
                                                     G4int At,  G4int Zt)
{
  const G4ParticleDefinition* hadron = aParticle->GetDefinition();
  G4double sumInelastic;
  G4int Nt = At - Zt;
  if(Nt < 0) Nt = 0;
  
  if( hadron == theKPlus )
  {
    sumInelastic =  GetHNinelasticXscVU(aParticle, At, Zt);
  }
  else
  {
    //sumInelastic  = Zt*GetHadronNucleonXscMK(aParticle, theProton);
    // sumInelastic += Nt*GetHadronNucleonXscMK(aParticle, theNeutron);    
    sumInelastic  = G4double(Zt)*GetHadronNucleonXscNS(aParticle, 1, 1);
    sumInelastic += G4double(Nt)*GetHadronNucleonXscNS(aParticle, 1, 0);    
  } 
  return sumInelastic;
}


/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon inelastic cross-section based on FTF-parametrisation 

G4double 
G4ComponentGGHadronNucleusXsc::GetHNinelasticXscVU(const G4DynamicParticle* aParticle, 
                                                 G4int At, G4int Zt)
{
  G4int PDGcode    = aParticle->GetDefinition()->GetPDGEncoding();
  G4int absPDGcode = std::abs(PDGcode);

  G4double Elab = aParticle->GetTotalEnergy();              
                          // (s - 2*0.88*GeV*GeV)/(2*0.939*GeV)/GeV;
  G4double Plab = aParticle->GetMomentum().mag();            
                          // std::sqrt(Elab * Elab - 0.88);

  Elab /= GeV;
  Plab /= GeV;

  G4double LogPlab    = G4Log( Plab );
  G4double sqrLogPlab = LogPlab * LogPlab;

  //G4cout<<"Plab = "<<Plab<<G4endl;

  G4double NumberOfTargetProtons = G4double(Zt); 
  G4double NumberOfTargetNucleons = G4double(At);
  G4double NumberOfTargetNeutrons = NumberOfTargetNucleons - NumberOfTargetProtons;

  if(NumberOfTargetNeutrons < 0.0) NumberOfTargetNeutrons = 0.0;

  G4double Xtotal, Xelastic, Xinelastic;

  if( absPDGcode > 1000 )  //------Projectile is baryon --------
  {
       G4double XtotPP = 48.0 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                         0.522*sqrLogPlab - 4.51*LogPlab;

       G4double XtotPN = 47.3 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                         0.513*sqrLogPlab - 4.27*LogPlab;

       G4double XelPP  = 11.9 + 26.9*G4Pow::GetInstance()->powA(Plab,-1.21) +
                         0.169*sqrLogPlab - 1.85*LogPlab;

       G4double XelPN  = 11.9 + 26.9*G4Pow::GetInstance()->powA(Plab,-1.21) +
                         0.169*sqrLogPlab - 1.85*LogPlab;

       Xtotal          = (NumberOfTargetProtons * XtotPP +
                          NumberOfTargetNeutrons * XtotPN);

       Xelastic        = (NumberOfTargetProtons * XelPP +
                          NumberOfTargetNeutrons * XelPN);
  }
  else if( PDGcode ==  211 ) //------Projectile is PionPlus -------
  {
       G4double XtotPiP = 16.4 + 19.3 *G4Pow::GetInstance()->powA(Plab,-0.42) +
                          0.19 *sqrLogPlab - 0.0 *LogPlab;

       G4double XtotPiN = 33.0 + 14.0 *G4Pow::GetInstance()->powA(Plab,-1.36) +
                          0.456*sqrLogPlab - 4.03*LogPlab;

       G4double XelPiP  =  0.0 + 11.4*G4Pow::GetInstance()->powA(Plab,-0.40) +
                           0.079*sqrLogPlab - 0.0 *LogPlab;

       G4double XelPiN  = 1.76 + 11.2*G4Pow::GetInstance()->powA(Plab,-0.64) +
                          0.043*sqrLogPlab - 0.0 *LogPlab;

       Xtotal           = ( NumberOfTargetProtons  * XtotPiP +
                            NumberOfTargetNeutrons * XtotPiN  );

       Xelastic         = ( NumberOfTargetProtons  * XelPiP  +
                            NumberOfTargetNeutrons * XelPiN   );
  }
  else if( PDGcode == -211 ) //------Projectile is PionMinus -------
  {
       G4double XtotPiP = 33.0 + 14.0 *G4Pow::GetInstance()->powA(Plab,-1.36) +
                          0.456*sqrLogPlab - 4.03*LogPlab;

       G4double XtotPiN = 16.4 + 19.3 *G4Pow::GetInstance()->powA(Plab,-0.42) +
                          0.19 *sqrLogPlab - 0.0 *LogPlab;

       G4double XelPiP  = 1.76 + 11.2*G4Pow::GetInstance()->powA(Plab,-0.64) +
                          0.043*sqrLogPlab - 0.0 *LogPlab;

       G4double XelPiN  =  0.0 + 11.4*G4Pow::GetInstance()->powA(Plab,-0.40) +
                           0.079*sqrLogPlab - 0.0 *LogPlab;

       Xtotal           = ( NumberOfTargetProtons  * XtotPiP +
                            NumberOfTargetNeutrons * XtotPiN  );

       Xelastic         = ( NumberOfTargetProtons  * XelPiP  +
                            NumberOfTargetNeutrons * XelPiN   );
  }
  else if( PDGcode ==  111 )  //------Projectile is PionZero  -------
  {
       G4double XtotPiP =(16.4 + 19.3 *G4Pow::GetInstance()->powA(Plab,-0.42) +
                          0.19 *sqrLogPlab - 0.0 *LogPlab +   //Pi+
                          33.0 + 14.0 *G4Pow::GetInstance()->powA(Plab,-1.36) +
                          0.456*sqrLogPlab - 4.03*LogPlab)/2; //Pi-

       G4double XtotPiN =(33.0 + 14.0 *G4Pow::GetInstance()->powA(Plab,-1.36) +
                          0.456*sqrLogPlab - 4.03*LogPlab +   //Pi+
                          16.4 + 19.3 *G4Pow::GetInstance()->powA(Plab,-0.42) +
                          0.19 *sqrLogPlab - 0.0 *LogPlab)/2; //Pi-

       G4double XelPiP  =( 0.0 + 11.4*G4Pow::GetInstance()->powA(Plab,-0.40) +
                           0.079*sqrLogPlab - 0.0 *LogPlab +    //Pi+
                           1.76 + 11.2*G4Pow::GetInstance()->powA(Plab,-0.64) +
                           0.043*sqrLogPlab - 0.0 *LogPlab)/2; //Pi-

       G4double XelPiN  =( 1.76 + 11.2*G4Pow::GetInstance()->powA(Plab,-0.64) +
                           0.043*sqrLogPlab - 0.0 *LogPlab +   //Pi+
                           0.0  + 11.4*G4Pow::GetInstance()->powA(Plab,-0.40) +
                           0.079*sqrLogPlab - 0.0 *LogPlab)/2; //Pi-

       Xtotal           = ( NumberOfTargetProtons  * XtotPiP +
                            NumberOfTargetNeutrons * XtotPiN  );

       Xelastic         = ( NumberOfTargetProtons  * XelPiP  +
                            NumberOfTargetNeutrons * XelPiN   );
  }
  else if( PDGcode == 321 ) //------Projectile is KaonPlus -------
  {
       G4double XtotKP = 18.1 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                         0.26 *sqrLogPlab - 1.0 *LogPlab;
       G4double XtotKN = 18.7 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                         0.21 *sqrLogPlab - 0.89*LogPlab;

       G4double XelKP  =  5.0 +  8.1*G4Pow::GetInstance()->powA(Plab,-1.8 ) +
                          0.16 *sqrLogPlab - 1.3 *LogPlab;

       G4double XelKN  =  7.3 +  0. *G4Pow::GetInstance()->powA(Plab,-0.  ) +
                          0.29 *sqrLogPlab - 2.4 *LogPlab;

       Xtotal          = ( NumberOfTargetProtons  * XtotKP +
                           NumberOfTargetNeutrons * XtotKN  );

       Xelastic        = ( NumberOfTargetProtons  * XelKP  +
                           NumberOfTargetNeutrons * XelKN   );
  }
  else if( PDGcode ==-321 )  //------Projectile is KaonMinus ------
  {
       G4double XtotKP = 32.1 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                         0.66 *sqrLogPlab - 5.6 *LogPlab;
       G4double XtotKN = 25.2 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                         0.38 *sqrLogPlab - 2.9 *LogPlab;

       G4double XelKP  =  7.3 +  0. *G4Pow::GetInstance()->powA(Plab,-0.  ) +
                          0.29 *sqrLogPlab - 2.4 *LogPlab;

       G4double XelKN  =  5.0 +  8.1*G4Pow::GetInstance()->powA(Plab,-1.8 ) +
                          0.16 *sqrLogPlab - 1.3 *LogPlab;

       Xtotal          = ( NumberOfTargetProtons  * XtotKP +
                           NumberOfTargetNeutrons * XtotKN  );

       Xelastic        = ( NumberOfTargetProtons  * XelKP  +
                           NumberOfTargetNeutrons * XelKN   );
  }
  else if( PDGcode == 311 ) //------Projectile is KaonZero ------
  {
       G4double XtotKP = ( 18.1 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                          0.26 *sqrLogPlab - 1.0 *LogPlab +   //K+
                          32.1 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                          0.66 *sqrLogPlab - 5.6 *LogPlab)/2; //K-

       G4double XtotKN = ( 18.7 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                          0.21 *sqrLogPlab - 0.89*LogPlab +   //K+
                          25.2 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                          0.38 *sqrLogPlab - 2.9 *LogPlab)/2; //K-

       G4double XelKP  = (  5.0 +  8.1*G4Pow::GetInstance()->powA(Plab,-1.8 )
                           + 0.16 *sqrLogPlab - 1.3 *LogPlab +   //K+
                           7.3 +  0. *G4Pow::GetInstance()->powA(Plab,-0.  ) +
                           0.29 *sqrLogPlab - 2.4 *LogPlab)/2; //K-

       G4double XelKN  = (  7.3 +  0. *G4Pow::GetInstance()->powA(Plab,-0.  ) +
                           0.29 *sqrLogPlab - 2.4 *LogPlab +   //K+
                           5.0 +  8.1*G4Pow::GetInstance()->powA(Plab,-1.8 ) +
                           0.16 *sqrLogPlab - 1.3 *LogPlab)/2; //K-

       Xtotal          = ( NumberOfTargetProtons  * XtotKP +
                           NumberOfTargetNeutrons * XtotKN  );

       Xelastic        = ( NumberOfTargetProtons  * XelKP  +
                           NumberOfTargetNeutrons * XelKN   );
  }
  else  //------Projectile is undefined, Nucleon assumed
  {
       G4double XtotPP = 48.0 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                         0.522*sqrLogPlab - 4.51*LogPlab;

       G4double XtotPN = 47.3 +  0. *G4Pow::GetInstance()->powA(Plab, 0.  ) +
                         0.513*sqrLogPlab - 4.27*LogPlab;

       G4double XelPP  = 11.9 + 26.9*G4Pow::GetInstance()->powA(Plab,-1.21) +
                         0.169*sqrLogPlab - 1.85*LogPlab;
       G4double XelPN  = 11.9 + 26.9*G4Pow::GetInstance()->powA(Plab,-1.21) +
                         0.169*sqrLogPlab - 1.85*LogPlab;

       Xtotal          = ( NumberOfTargetProtons  * XtotPP +
                           NumberOfTargetNeutrons * XtotPN  );

       Xelastic        = ( NumberOfTargetProtons  * XelPP  +
                           NumberOfTargetNeutrons * XelPN   );
  }
  Xinelastic = Xtotal - Xelastic;

  if( Xinelastic < 0.) Xinelastic = 0.;

  return Xinelastic*= millibarn;
}

////////////////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4ComponentGGHadronNucleusXsc::GetNucleusRadius(const G4DynamicParticle* , 
                                              const G4Element* anElement)
{
  G4int At = G4lrint(anElement->GetN());
  G4double oneThird = 1.0/3.0;
  G4double cubicrAt = G4Pow::GetInstance()->powA(G4double(At), oneThird); 

  G4double R;  // = fRadiusConst*cubicrAt;
  /*  
  G4double tmp = G4Pow::GetInstance()->powA( cubicrAt-1., 3.);
  tmp         += At;
  tmp         *= 0.5;

  if (At > 20.)   // 20.
  {
    R = fRadiusConst*G4Pow::GetInstance()->powA (tmp, oneThird); 
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

  if (At > 20)   // 20.
  {
    R *= ( a1 + b1*G4Exp( -(At - meanA)/tauA1) ); 
  }
  else if (At > 3)
  {
    R *= ( 1.0 + b2*( 1. - G4Exp( (At - meanA)/tauA2) ) ); 
  }
  else 
  {
    R *= ( 1.0 + b3*( 1. - G4Exp( (At - meanA)/tauA3) ) ); 
  }  
  return R;
 
}
////////////////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4ComponentGGHadronNucleusXsc::GetNucleusRadius(G4int At)
{
  G4double oneThird = 1.0/3.0;
  G4double cubicrAt = G4Pow::GetInstance()->powA(G4double(At), oneThird); 

  G4double R;  // = fRadiusConst*cubicrAt;

  /*
  G4double tmp = G4Pow::GetInstance()->powA( cubicrAt-1., 3.);
  tmp         += At;
  tmp         *= 0.5;

  if (At > 20.)
  {
    R = fRadiusConst*G4Pow::GetInstance()->powA (tmp, oneThird); 
  }
  else
  {
    R = fRadiusConst*cubicrAt; 
  }
  */

  R = fRadiusConst*cubicrAt;

  G4double meanA = 20.;
  G4double tauA  = 20.; 

  if (At > 20)   // 20.
  {
    R *= ( 0.8 + 0.2*G4Exp( -(G4double(At) - meanA)/tauA) ); 
  }
  else
  {
    R *= ( 1.0 + 0.1*( 1. - G4Exp( (G4double(At) - meanA)/tauA) ) ); 
  }

  return R;
}

////////////////////////////////////////////////////////////////////////////////////
//
//

G4double G4ComponentGGHadronNucleusXsc::CalculateEcmValue( const G4double mp , 
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

G4double G4ComponentGGHadronNucleusXsc::CalcMandelstamS( const G4double mp , 
                                                       const G4double mt , 
                                                       const G4double Plab )
{
  G4double Elab = std::sqrt ( mp * mp + Plab * Plab );
  G4double sMand  = mp*mp + mt*mt + 2*Elab*mt ;

  return sMand;
}

////////////////////////////////////////////////////////////////////////////////////
//
//

void G4ComponentGGHadronNucleusXsc::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4ComponentGGHadronNucleusXsc calculates total, inelastic and\n"
          << "elastic cross sections for hadron-nucleus cross sections using\n"
          << "the Glauber model with Gribov corrections.  It is valid for all\n"
          << "targets except hydrogen, and for incident p, pbar, n, sigma-,\n"
          << "pi+, pi-, K+, K- and gammas with energies above 3 GeV.  This is\n"
          << "a cross section component which is to be used to build a cross\n"
          << "data set.\n";
}


///////////////////////////////////////////////////////////////////////////////
//
// Correction arrays for GG <-> Bar changea at ~ 90 GeV

const G4double G4ComponentGGHadronNucleusXsc::fNeutronBarCorrectionTot[93] = {

  1.0, 1.0,     1.42517e+00,  // 1.118517e+00, 
1.082002e+00, 1.116171e+00, 1.078747e+00, 1.061315e+00, 
1.058205e+00, 1.082663e+00, 1.068500e+00, 1.076912e+00, 1.083475e+00, 1.079117e+00, 
1.071856e+00, 1.071990e+00, 1.073774e+00, 1.079356e+00, 1.081314e+00, 1.082056e+00,
1.090772e+00, 1.096776e+00, 1.095828e+00, 1.097678e+00, 1.099157e+00, 1.103677e+00, 
1.105132e+00, 1.109806e+00, 1.110816e+00, 1.117378e+00, 1.115165e+00, 1.115710e+00, 
1.111855e+00, 1.110482e+00, 1.110112e+00, 1.106676e+00, 1.108706e+00, 1.105549e+00, 
1.106318e+00, 1.106242e+00, 1.107672e+00, 1.107342e+00, 1.108119e+00, 1.106655e+00, 
1.102588e+00, 1.096657e+00, 1.092920e+00, 1.086629e+00, 1.083592e+00, 1.076030e+00, 
1.083777e+00, 1.089460e+00, 1.086545e+00, 1.079924e+00, 1.082218e+00, 1.077798e+00, 
1.077062e+00, 1.072825e+00, 1.072241e+00, 1.072104e+00, 1.072490e+00, 1.069829e+00, 
1.070398e+00, 1.065458e+00, 1.064968e+00, 1.060524e+00, 1.060048e+00, 1.057620e+00, 
1.056428e+00, 1.055366e+00, 1.055017e+00, 1.052304e+00, 1.051767e+00, 1.049728e+00, 
1.048745e+00, 1.047399e+00, 1.045876e+00, 1.042972e+00, 1.041824e+00, 1.039993e+00, 
1.039021e+00, 1.036627e+00, 1.034176e+00, 1.032526e+00, 1.033633e+00, 1.036107e+00, 
1.037803e+00, 1.031266e+00, 1.032991e+00, 1.033284e+00, 1.035015e+00, 1.033945e+00, 
1.037075e+00, 1.034721e+00

};

const G4double G4ComponentGGHadronNucleusXsc::fNeutronBarCorrectionIn[93] = {

1.0, 1.0,     1.167421e+00, 1.156250e+00, 1.205364e+00, 1.154225e+00, 1.120391e+00, // 6
1.124632e+00, 1.129460e+00, 1.107863e+00, 1.102152e+00, 1.104593e+00, 1.100285e+00, // 12
1.098450e+00, 1.092677e+00, 1.101124e+00, 1.106461e+00, 1.115049e+00, 1.123903e+00, // 18
1.126661e+00, 1.131259e+00, 1.133949e+00, 1.134185e+00, 1.133767e+00, 1.132813e+00, // 24
1.131515e+00, 1.144338e+00, // 1.130338e+00, 
1.134171e+00, 1.139206e+00, 1.148474e+00, // 1.141474e+00, 
1.142189e+00, 
1.140725e+00, 1.140100e+00, 1.139848e+00, 1.137674e+00, 1.138645e+00, 1.136339e+00, 
1.136439e+00, 1.135946e+00, 1.136431e+00, 1.135702e+00, 1.135703e+00, 1.134113e+00, 
1.131935e+00, 1.128381e+00, 1.126373e+00, 1.122453e+00, 1.120908e+00, 1.115953e+00, 
1.115947e+00, 1.114426e+00, 1.111749e+00, 1.106207e+00, 1.107494e+00, 1.103622e+00, 
1.102576e+00, 1.098816e+00, 1.097889e+00, 1.097306e+00, 1.097130e+00, 1.094578e+00, 
1.094552e+00, 1.090222e+00, 1.089358e+00, 1.085409e+00, 1.084560e+00, 1.082182e+00, 
1.080773e+00, 1.079464e+00, 1.078724e+00, 1.076121e+00, 1.075235e+00, 1.073159e+00, 
1.071920e+00, 1.070395e+00, 1.069503e+00, 1.067525e+00, 1.066919e+00, 1.065779e+00, 
1.065319e+00, 1.063730e+00, 1.062092e+00, 1.061085e+00, 1.059908e+00, 1.059815e+00, 
1.059109e+00, 1.051920e+00, 1.051258e+00, 1.049473e+00, 1.048823e+00, 1.045984e+00, 
1.046435e+00, 1.042614e+00

};

const G4double G4ComponentGGHadronNucleusXsc::fProtonBarCorrectionTot[93] = {

1.0, 1.0,     
1.118515e+00, 1.082000e+00, 1.116169e+00, 1.078745e+00, 1.061313e+00, 1.058203e+00, 
1.082661e+00, 1.068498e+00, 1.076910e+00, 1.083474e+00, 1.079115e+00, 1.071854e+00, 
1.071988e+00, 1.073772e+00, 1.079355e+00, 1.081312e+00, 1.082054e+00, 1.090770e+00, 
1.096774e+00, 1.095827e+00, 1.097677e+00, 1.099156e+00, 1.103676e+00, 1.105130e+00, 
1.109805e+00, 1.110814e+00, 1.117377e+00, 1.115163e+00, 1.115708e+00, 1.111853e+00, 
1.110480e+00, 1.110111e+00, 1.106674e+00, 1.108705e+00, 1.105548e+00, 1.106317e+00, 
1.106241e+00, 1.107671e+00, 1.107341e+00, 1.108118e+00, 1.106654e+00, 1.102586e+00, 
1.096655e+00, 1.092918e+00, 1.086628e+00, 1.083590e+00, 1.076028e+00, 1.083776e+00, 
1.089458e+00, 1.086543e+00, 1.079923e+00, 1.082216e+00, 1.077797e+00, 1.077061e+00, 
1.072824e+00, 1.072239e+00, 1.072103e+00, 1.072488e+00, 1.069828e+00, 1.070396e+00, 
1.065456e+00, 1.064966e+00, 1.060523e+00, 1.060047e+00, 1.057618e+00, 1.056427e+00, 
1.055365e+00, 1.055016e+00, 1.052303e+00, 1.051766e+00, 1.049727e+00, 1.048743e+00, 
1.047397e+00, 1.045875e+00, 1.042971e+00, 1.041823e+00, 1.039992e+00, 1.039019e+00, 
1.036626e+00, 1.034175e+00, 1.032525e+00, 1.033632e+00, 1.036106e+00, 1.037802e+00, 
1.031265e+00, 1.032990e+00, 1.033283e+00, 1.035014e+00, 1.033944e+00, 1.037074e+00, 
1.034720e+00 

};

const G4double G4ComponentGGHadronNucleusXsc::fProtonBarCorrectionIn[93] = {

1.0, 1.0,     
1.147419e+00, // 1.167419e+00, 
1.156248e+00, 1.205362e+00, 1.154224e+00, 1.120390e+00, 1.124630e+00, // 7 
1.129459e+00, 1.107861e+00, 1.102151e+00, 1.104591e+00, 1.100284e+00, 1.098449e+00, // 13
1.092675e+00, 1.101122e+00, 1.106460e+00, 1.115048e+00, 1.123902e+00, 1.126659e+00, // 19
1.131258e+00, 1.133948e+00, 1.134183e+00, 1.133766e+00, 1.132812e+00, 1.131514e+00, // 25
1.140337e+00, // 1.130337e+00, 

1.134170e+00, 1.139205e+00, 1.151472e+00,  // 1.141472e+00, 
1.142188e+00, 1.140724e+00, 
1.140099e+00, 1.139847e+00, 1.137672e+00, 1.138644e+00, 1.136338e+00, 1.136438e+00, 
1.135945e+00, 1.136429e+00, 1.135701e+00, 1.135702e+00, 1.134112e+00, 1.131934e+00, 
1.128380e+00, 1.126371e+00, 1.122452e+00, 1.120907e+00, 1.115952e+00, 1.115946e+00, 
1.114425e+00, 1.111748e+00, 1.106205e+00, 1.107493e+00, 1.103621e+00, 1.102575e+00, 
1.098815e+00, 1.097888e+00, 1.097305e+00, 1.097129e+00, 1.094577e+00, 1.094551e+00, 
1.090221e+00, 1.089357e+00, 1.085408e+00, 1.084559e+00, 1.082181e+00, 1.080772e+00, 
1.079463e+00, 1.078723e+00, 1.076120e+00, 1.075234e+00, 1.073158e+00, 1.071919e+00, 
1.070394e+00, 1.069502e+00, 1.067524e+00, 1.066918e+00, 1.065778e+00, 1.065318e+00, 
1.063729e+00, 1.062091e+00, 1.061084e+00, 1.059907e+00, 1.059814e+00, 1.059108e+00, 
1.051919e+00, 1.051257e+00, 1.049472e+00, 1.048822e+00, 1.045983e+00, 1.046434e+00, 
1.042613e+00 

};


const G4double G4ComponentGGHadronNucleusXsc::fPionPlusBarCorrectionTot[93] = {

1.0, 1.0,     
1.075927e+00, 1.074407e+00, 1.126098e+00, 1.100127e+00, 1.089742e+00, 1.083536e+00, 
1.089988e+00, 1.103566e+00, 1.096922e+00, 1.126573e+00, 1.132734e+00, 1.136512e+00, 
1.136629e+00, 1.133086e+00, 1.132428e+00, 1.129299e+00, 1.125622e+00, 1.126992e+00, 
1.127840e+00, 1.162670e+00, 1.160392e+00, 1.157864e+00, 1.157227e+00, 1.154627e+00, 
1.192555e+00, 1.197243e+00, 1.197911e+00, 1.200326e+00, 1.220053e+00, 1.215019e+00, 
1.211703e+00, 1.209080e+00, 1.204248e+00, 1.203328e+00, 1.198671e+00, 1.196840e+00, 
1.194392e+00, 1.193037e+00, 1.190408e+00, 1.188583e+00, 1.206127e+00, 1.210028e+00, 
1.206434e+00, 1.204456e+00, 1.200547e+00, 1.199058e+00, 1.200174e+00, 1.200276e+00, 
1.198912e+00, 1.213048e+00, 1.207160e+00, 1.208020e+00, 1.203814e+00, 1.202380e+00, 
1.198306e+00, 1.197002e+00, 1.196027e+00, 1.195449e+00, 1.192563e+00, 1.192135e+00, 
1.187556e+00, 1.186308e+00, 1.182124e+00, 1.180900e+00, 1.178224e+00, 1.176471e+00, 
1.174811e+00, 1.173702e+00, 1.170827e+00, 1.169581e+00, 1.167205e+00, 1.165626e+00, 
1.180244e+00, 1.177626e+00, 1.175121e+00, 1.173903e+00, 1.172192e+00, 1.171128e+00, 
1.168997e+00, 1.166826e+00, 1.164130e+00, 1.165412e+00, 1.165504e+00, 1.165020e+00, 
1.158462e+00, 1.158014e+00, 1.156519e+00, 1.156081e+00, 1.153602e+00, 1.154190e+00, 
1.152974e+00
 
};

const G4double G4ComponentGGHadronNucleusXsc::fPionPlusBarCorrectionIn[93] = {

1.0, 1.0,    
1.140246e+00, 1.097872e+00, 1.104301e+00, 1.068722e+00, 1.056495e+00, 1.062622e+00, // 7
1.047987e+00, 1.037032e+00, 1.035686e+00, 1.042870e+00, 1.052222e+00, 1.075100e+00, // 13
1.084480e+00, 1.078286e+00, 1.081488e+00, 1.089713e+00, 1.099105e+00, 1.098003e+00, // 19
1.102175e+00, 1.117707e+00, 1.121734e+00, 1.125229e+00, 1.126457e+00, 1.128905e+00, // 25
1.163312e+00, 1.126263e+00, 1.126459e+00, 1.135191e+00, 1.116986e+00, 1.117184e+00, // 31
1.117037e+00, 1.116777e+00, 1.115858e+00, 1.115745e+00, 1.114489e+00, 1.113993e+00, // 37
1.113226e+00, 1.112818e+00, 1.111890e+00, 1.111238e+00, 1.111209e+00, 1.111775e+00, // 43
1.110256e+00, 1.109414e+00, 1.107647e+00, 1.106980e+00, 1.106096e+00, 1.107331e+00, // 49
1.107849e+00, 1.106407e+00, 1.103426e+00, 1.103896e+00, 1.101756e+00, 1.101031e+00, // 55
1.098915e+00, 1.098260e+00, 1.097768e+00, 1.097487e+00, 1.095964e+00, 1.095773e+00, // 61
1.093348e+00, 1.092687e+00, 1.090465e+00, 1.089821e+00, 1.088394e+00, 1.087462e+00, // 67
1.086571e+00, 1.085997e+00, 1.084451e+00, 1.083798e+00, 1.082513e+00, 1.081670e+00, // 73
1.080735e+00, 1.075659e+00, 1.074341e+00, 1.073689e+00, 1.072787e+00, 1.072237e+00, // 79
1.071107e+00, 1.069955e+00, 1.074856e+00, 1.065873e+00, 1.065938e+00, 1.065694e+00, 
1.062192e+00, 1.061967e+00, 1.061180e+00, 1.060960e+00, 1.059646e+00, 1.059975e+00, 
1.059658e+00
 
};


const G4double G4ComponentGGHadronNucleusXsc::fPionMinusBarCorrectionTot[93] = {

1.0, 1.0,     
1.3956e+00, 1.077959e+00, 1.129145e+00, 1.102088e+00, 1.089765e+00, 1.083542e+00,  // 7
1.089995e+00, 1.104895e+00, 1.097154e+00, 1.127663e+00, 1.133063e+00, 1.137425e+00, // 13
1.136724e+00, 1.133859e+00, 1.132498e+00, 1.130276e+00, 1.127896e+00, 1.127656e+00, // 19
1.127905e+00, 1.164210e+00, 1.162259e+00, 1.160075e+00, 1.158978e+00, 1.156649e+00, // 25 
1.194157e+00, 1.199177e+00, 1.198983e+00, 1.202325e+00, 1.221967e+00, 1.217548e+00, 
1.214389e+00, 1.211760e+00, 1.207335e+00, 1.206081e+00, 1.201766e+00, 1.199779e+00, 
1.197283e+00, 1.195706e+00, 1.193071e+00, 1.191115e+00, 1.208838e+00, 1.212681e+00, 
1.209235e+00, 1.207163e+00, 1.203451e+00, 1.201807e+00, 1.203283e+00, 1.203388e+00, 
1.202244e+00, 1.216509e+00, 1.211066e+00, 1.211504e+00, 1.207539e+00, 1.205991e+00, 
1.202143e+00, 1.200724e+00, 1.199595e+00, 1.198815e+00, 1.196025e+00, 1.195390e+00, 
1.191137e+00, 1.189791e+00, 1.185888e+00, 1.184575e+00, 1.181996e+00, 1.180229e+00, 
1.178545e+00, 1.177355e+00, 1.174616e+00, 1.173312e+00, 1.171016e+00, 1.169424e+00, 
1.184120e+00, 1.181478e+00, 1.179085e+00, 1.177817e+00, 1.176124e+00, 1.175003e+00, 
1.172947e+00, 1.170858e+00, 1.168170e+00, 1.169397e+00, 1.169304e+00, 1.168706e+00, 
1.162774e+00, 1.162217e+00, 1.160740e+00, 1.160196e+00, 1.157857e+00, 1.158220e+00, 
1.157267e+00 
};


const G4double G4ComponentGGHadronNucleusXsc::fPionMinusBarCorrectionIn[93] = {

1.0, 1.0,    
1.463e+00,    1.100898e+00, 1.106773e+00, 1.070289e+00, 1.040514e+00, 1.062628e+00, // 7
1.047992e+00, 1.038041e+00, 1.035862e+00, 1.043679e+00, 1.052466e+00, 1.065780e+00, // 13
1.070551e+00, 1.078869e+00, 1.081541e+00, 1.090455e+00, 1.100847e+00, 1.098511e+00, // 19 
1.102226e+00, 1.118865e+00, 1.123143e+00, 1.126904e+00, 1.127785e+00, 1.130444e+00, // 25
1.148502e+00, 1.127678e+00, 1.127244e+00, 1.123634e+00, 1.118347e+00, 1.118988e+00, 
1.118957e+00, 1.118696e+00, 1.118074e+00, 1.117722e+00, 1.116717e+00, 1.116111e+00, 
1.115311e+00, 1.114745e+00, 1.113814e+00, 1.113069e+00, 1.113141e+00, 1.113660e+00, 
1.112249e+00, 1.111343e+00, 1.109718e+00, 1.108942e+00, 1.108310e+00, 1.109549e+00, 
1.110227e+00, 1.108846e+00, 1.106183e+00, 1.106354e+00, 1.104388e+00, 1.103583e+00, 
1.101632e+00, 1.100896e+00, 1.100296e+00, 1.099873e+00, 1.098420e+00, 1.098082e+00, 
1.095892e+00, 1.095162e+00, 1.093144e+00, 1.092438e+00, 1.091083e+00, 1.090142e+00, 
1.089236e+00, 1.088604e+00, 1.087159e+00, 1.086465e+00, 1.085239e+00, 1.084388e+00, 
1.083473e+00, 1.078373e+00, 1.077136e+00, 1.076450e+00, 1.075561e+00, 1.074973e+00, 
1.073898e+00, 1.072806e+00, 1.067706e+00, 1.068684e+00, 1.068618e+00, 1.068294e+00, 
1.065241e+00, 1.064939e+00, 1.064166e+00, 1.063872e+00, 1.062659e+00, 1.062828e+00, 
1.062699e+00 

};


//
//
///////////////////////////////////////////////////////////////////////////////////////
