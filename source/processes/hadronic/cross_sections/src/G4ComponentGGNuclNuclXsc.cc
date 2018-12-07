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
// 24.11.08 V. Grichine - first implementation
//
// 04.09.18 V. Ivantchenko Major revision of interfaces and implementation

#include "G4ComponentGGNuclNuclXsc.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadTmpUtil.hh"
#include "G4HadronNucleonXsc.hh"
#include "G4ComponentGGHadronNucleusXsc.hh" 
#include "G4Pow.hh"

static const G4double inve = 1./CLHEP::eplus;

G4ComponentGGNuclNuclXsc::G4ComponentGGNuclNuclXsc() 
 : G4VComponentCrossSection("Glauber-Gribov Nucl-nucl"),
   fRadiusConst(1.08*fermi),  // 1.1, 1.3 ?
   fTotalXsc(0.0), fElasticXsc(0.0), fInelasticXsc(0.0), fProductionXsc(0.0),
   fDiffractionXsc(0.0), fParticle(nullptr), fEnergy(0.0), fZ(0), fA(0)
{
  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  fHNXsc = new G4HadronNucleonXsc();
  fHadrNucl = new G4ComponentGGHadronNucleusXsc(); 
  fNist  = G4NistManager::Instance();
  fCalc  = G4Pow::GetInstance();
}

G4ComponentGGNuclNuclXsc::~G4ComponentGGNuclNuclXsc()
{
  delete fHNXsc;
}

//////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetTotalElementCrossSection(
         const G4ParticleDefinition* aParticle, G4double kinEnergy, 
	 G4int Z, G4double A)
{
  ComputeCrossSections(aParticle, kinEnergy, Z, G4lrint(A));
  return fTotalXsc;
}

////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetTotalIsotopeCrossSection(
         const G4ParticleDefinition* aParticle, G4double kinEnergy,
	 G4int Z, G4int A)
{
  ComputeCrossSections(aParticle, kinEnergy, Z, A);
  return fTotalXsc;
}

/////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetInelasticElementCrossSection(
         const G4ParticleDefinition* aParticle, G4double kinEnergy, 
	 G4int Z, G4double A)
{
  ComputeCrossSections(aParticle, kinEnergy, Z, G4lrint(A));
  return fInelasticXsc;
}

////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetInelasticIsotopeCrossSection(
         const G4ParticleDefinition* aParticle, G4double kinEnergy, 
	 G4int Z, G4int A)
{
  ComputeCrossSections(aParticle, kinEnergy, Z, A);
  return fInelasticXsc;
}

//////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetElasticElementCrossSection(
         const G4ParticleDefinition* aParticle, G4double kinEnergy, 
	 G4int Z, G4double A)
{
  ComputeCrossSections(aParticle, kinEnergy, Z, G4lrint(A));
  return fElasticXsc;
}

///////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetElasticIsotopeCrossSection(
         const G4ParticleDefinition* aParticle, G4double kinEnergy, 
	 G4int Z, G4int A)
{
  ComputeCrossSections(aParticle, kinEnergy, Z, A);
  return fElasticXsc;
}

////////////////////////////////////////////////////////////////
 
G4double G4ComponentGGNuclNuclXsc::ComputeQuasiElasticRatio(
         const G4ParticleDefinition* aParticle, G4double kinEnergy, 
	 G4int Z, G4int A)
{
  ComputeCrossSections(aParticle, kinEnergy, Z, A);
  return (fInelasticXsc > fProductionXsc) 
    ? (fInelasticXsc - fProductionXsc)/fInelasticXsc : 0.0;
}

//////////////////////////////////////////////////////////////////////
 
void G4ComponentGGNuclNuclXsc::BuildPhysicsTable(const G4ParticleDefinition&)
{}

//////////////////////////////////////////////////////////////////////

void G4ComponentGGNuclNuclXsc::DumpPhysicsTable(const G4ParticleDefinition&)
{
  G4cout << "G4ComponentGGNuclNuclXsc: uses Glauber-Gribov formula" << G4endl;
}

//////////////////////////////////////////////////////////////////////

void G4ComponentGGNuclNuclXsc::Description(std::ostream& outFile) const
{
  outFile << "G4ComponentGGNuclNuclXsc calculates total, inelastic and\n"
          << "elastic cross sections for nucleus-nucleus collisions using\n"
          << "the Glauber model with Gribov corrections.  It is valid for\n"
          << "all incident energies above 100 keV./n"
	  << "For the hydrogen target G4HadronNucleonXsc class is used.\n";
}

/////////////////////////////////////////////////////////////////////

G4bool G4ComponentGGNuclNuclXsc::IsElementApplicable(const G4DynamicParticle*, 
					             G4int, const G4Material*)
{
  return true;
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculates total and inelastic Xsc, derives elastic as total - inelastic 
// accordong to Glauber model with Gribov correction calculated in the dipole 
// approximation on light cone. Gaussian density of point-like nucleons helps 
// to calculate rest integrals of the model. [1] B.Z. Kopeliovich,
// nucl-th/0306044 + simplification above

void G4ComponentGGNuclNuclXsc::ComputeCrossSections(
     const G4ParticleDefinition* aParticle, G4double kinEnergy, 
     G4int Z, G4int A)
{
  // check cache
  if(aParticle == fParticle && fZ == Z && fA == A && kinEnergy == fEnergy)
    { return; }
  fParticle = aParticle;
  fZ = Z;
  fA = A;
  fEnergy = kinEnergy;

  G4int pZ = G4lrint(aParticle->GetPDGCharge()*inve);
  G4int pA = aParticle->GetBaryonNumber();

  // hydrogen
  if(1 == Z && 1 == A) {
    G4double e = kinEnergy*CLHEP::proton_mass_c2/aParticle->GetPDGMass();
    fHadrNucl->ComputeCrossSections(theProton, e, pZ, pA);
    fTotalXsc = fHadrNucl->GetTotalGlauberGribovXsc(); 
    fElasticXsc = fHadrNucl->GetElasticGlauberGribovXsc(); 
    fInelasticXsc = fHadrNucl->GetInelasticGlauberGribovXsc(); 
    fProductionXsc = fHadrNucl->GetProductionGlauberGribovXsc(); 
    fDiffractionXsc = fHadrNucl->GetDiffractionGlauberGribovXsc(); 
    return;
  }
  static const G4double cofInelastic = 2.4;
  static const G4double cofTotal = 2.0;

  G4double pTkin = kinEnergy/(G4double)pA;

  G4int pN = pA - pZ;
  G4int tN = A - Z;

  G4double tR = GetNucleusRadius(  Z,  A);  
  G4double pR = GetNucleusRadius( pZ, pA); 

  G4double cB = ComputeCoulombBarier(aParticle, kinEnergy, Z, A, pR, tR);

  if ( cB > 0. ) 
  {
    G4double sigma = (pZ*Z+pN*tN)*fHNXsc->HadronNucleonXscNS(theProton, theProton, pTkin);
    G4double ppInXsc = fHNXsc->GetInelasticHadronNucleonXsc();

    sigma += (pZ*tN+pN*Z)*fHNXsc->HadronNucleonXscNS(theNeutron, theProton, pTkin);
    G4double npInXsc = fHNXsc->GetInelasticHadronNucleonXsc();

    // G4cout<<"ppInXsc = "<<ppInXsc/millibarn<<"; npInXsc = "<<npInXsc/millibarn<<G4endl;
    // G4cout<<"npTotXsc = "<<fHNXsc->GetTotalHadronNucleonXsc()/millibarn<<"; npElXsc = "
    //                      <<fHNXsc->GetElasticHadronNucleonXsc()/millibarn<<G4endl;

    G4double nucleusSquare = cofTotal*CLHEP::pi*( pR*pR + tR*tR ); // basically 2piRR

    G4double ratio= sigma/nucleusSquare;
    fTotalXsc     = nucleusSquare*G4Log( 1. + ratio )*cB;
    fInelasticXsc = nucleusSquare*G4Log( 1. + cofInelastic*ratio )*cB/cofInelastic;
    fElasticXsc   = std::max(fTotalXsc - fInelasticXsc, 0.0);

    G4double difratio = ratio/(1.+ratio);
    fDiffractionXsc = 0.5*nucleusSquare*( difratio - G4Log( 1. + difratio ) );

    G4double xratio= ((pZ*Z+pN*tN)*ppInXsc + (pZ*tN+pN*Z)*npInXsc)/nucleusSquare;
    fProductionXsc = nucleusSquare*G4Log( 1. + cofInelastic*xratio)*cB/cofInelastic;
    fProductionXsc = std::min(fProductionXsc, fInelasticXsc);
  }
  else
  {
    fInelasticXsc  = 0.;
    fTotalXsc      = 0.;
    fElasticXsc    = 0.;
    fProductionXsc = 0.;
    fDiffractionXsc= 0.;
  }
}

///////////////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::ComputeCoulombBarier(
                     const G4ParticleDefinition* aParticle,
		     G4double pTkin, G4int Z, G4int A,
		     G4double pR, G4double tR)
{
  G4int pZ = aParticle->GetPDGCharge()*inve;
  G4double pM = aParticle->GetPDGMass();
  G4double tM = G4NucleiProperties::GetNuclearMass(A, Z); 
  G4double pElab = pTkin + pM;
  G4double totEcm = std::sqrt(pM*pM + tM*tM + 2.*pElab*tM);
  G4double totTcm = totEcm - pM -tM;

  static const G4double qfact = CLHEP::fine_structure_const*CLHEP::hbarc; 
  G4double bC = qfact*pZ*Z*0.5/(pR + tR);

  G4double ratio = (totTcm <= bC ) ? 0. : 1. - bC/totTcm;
  // G4cout<<"G4ComponentGGNuclNuclXsc::ComputeCoulombBarier= "<<ratio
  // <<"; pTkin(GeV)= " <<pTkin/GeV<<"; 
  // " pPlab = "<<pPlab/GeV<<"; bC = "<<bC/GeV<<"; pTcm = "
  // <<pTcm/GeV<<G4endl;
  return ratio;
}

//////////////////////////////////////////////////////////////////////////
//
// Return single-diffraction/inelastic cross-section ratio

G4double G4ComponentGGNuclNuclXsc::GetRatioSD(
         const G4DynamicParticle* aParticle, G4double tA, G4double tZ)
{
  ComputeCrossSections(aParticle->GetDefinition(), 
                       aParticle->GetKineticEnergy(), 
		       G4lrint(tZ), G4lrint(tA));

  return (fInelasticXsc > 0.0) ? fDiffractionXsc/fInelasticXsc : 0.0;
}

//////////////////////////////////////////////////////////////////////////
//
// Return quasi-elastic/inelastic cross-section ratio

G4double G4ComponentGGNuclNuclXsc::GetRatioQE(
         const G4DynamicParticle* aParticle, G4double tA, G4double tZ)
{
  ComputeCrossSections(aParticle->GetDefinition(), 
                       aParticle->GetKineticEnergy(), 
		       G4lrint(tZ), G4lrint(tA));

  return (fInelasticXsc > 0.0) ? 1.0 - fProductionXsc/fInelasticXsc : 0.0;
}

///////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to differnt parametrisations:
// [2] E. Levin, hep-ph/9710546
// [3] U. Dersch, et al, hep-ex/9910052
// [4] M.J. Longo, et al, Phys.Rev.Lett. 33 (1974) 725 

G4double 
G4ComponentGGNuclNuclXsc::GetHadronNucleonXsc(const G4DynamicParticle* aParticle, 
                                              const G4Element* anElement)
{
  G4int At = G4lrint(anElement->GetN());  // number of nucleons 
  G4int Zt = anElement->GetZasInt();      // number of protons
  return GetHadronNucleonXsc(aParticle, At, Zt);
}

///////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to differnt parametrisations:
// [2] E. Levin, hep-ph/9710546
// [3] U. Dersch, et al, hep-ex/9910052
// [4] M.J. Longo, et al, Phys.Rev.Lett. 33 (1974) 725 

G4double 
G4ComponentGGNuclNuclXsc::GetHadronNucleonXsc(const G4DynamicParticle* aParticle, 
                                              G4int At, G4int Zt)
{
  return fHadrNucl->GetHadronNucleonXsc(aParticle, At, Zt);
}

///////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to PDG parametrisation (2005):
// http://pdg.lbl.gov/2006/reviews/hadronicrpp.pdf
//  At = number of nucleons,  Zt = number of protons 

G4double 
G4ComponentGGNuclNuclXsc::GetHadronNucleonXscPDG(const G4ParticleDefinition* pParticle, 
                                                 G4double pTkin, 
                                                 const G4ParticleDefinition* tParticle)
{
  G4double res = 0.0;
  if(tParticle == theProton) {
    res = fHNXsc->HadronNucleonXscPDG(pParticle, theProton, pTkin); 
  } else if(tParticle == theNeutron) {
    res = fHNXsc->HadronNucleonXscPDG(pParticle, theNeutron, pTkin); 
  } else {
    G4int Zt = tParticle->GetAtomicNumber();
    G4int At = tParticle->GetAtomicMass();
    fHadrNucl->ComputeCrossSections(pParticle, pTkin, Zt, At);
    res = fHadrNucl->GetTotalGlauberGribovXsc();
  }
  return res;
}

///////////////////////////////////////////////////////////////////////////////
//
// Returns total nucleon-nucleon cross-section based on N. Starkov parametrisation 
// of data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database
// projectile nucleon is pParticle with pTkin shooting target nucleon tParticle

G4double 
G4ComponentGGNuclNuclXsc::GetHadronNucleonXscNS(const G4ParticleDefinition* pParticle, 
						G4double pTkin, 
						const G4ParticleDefinition* tParticle)
{
  G4int Zt = 1;
  G4int At = 1;
  if(tParticle == theNeutron) { Zt = 0; }
  else if(tParticle != theProton) {
    Zt = tParticle->GetAtomicNumber();
    At = tParticle->GetAtomicMass();
  }
  fHadrNucl->ComputeCrossSections(pParticle, pTkin, Zt, At);
  return fHadrNucl->GetTotalGlauberGribovXsc();
}

/////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon inelastic cross-section based on FTF-parametrisation 

G4double 
G4ComponentGGNuclNuclXsc::GetHNinelasticXscVU(const G4DynamicParticle* aParticle, 
                                              G4int At, G4int Zt)
{
  return fHadrNucl->GetHNinelasticXscVU(aParticle, At, Zt);
}

///////////////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetNucleusRadius(const G4DynamicParticle*, 
                                                    const G4Element* anElement)
{
  G4double At = anElement->GetN();
  G4double R = fRadiusConst*fCalc->A13(At);

  static const G4double meanA  = 21.;
  static const G4double tauA1  = 40.; 
  static const G4double tauA2  = 10.; 
  static const G4double tauA3  = 5.; 

  static const G4double a1 = 0.85;
  static const G4double b1 = 1. - a1;

  static const G4double b2 = 0.3;
  static const G4double b3 = 4.;

  if (At > 20.)   // 20.
  {
    R *= ( a1 + b1*G4Exp( -(At - meanA)/tauA1) ); 
  }
  else if (At > 3.5)
  {
    R *= ( 1.0 + b2*( 1. - G4Exp( (At - meanA)/tauA2) ) ); 
  }
  else 
  {
    R *= ( 1.0 + b3*( 1. - G4Exp( (At - meanA)/tauA3) ) ); 
  }
  return R;
}

///////////////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetNucleusRadius(G4int Zt, G4int At)
{
  return GetNucleusRadiusDE(Zt, At);
}

///////////////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetNucleusRadiusGG(G4int At)
{
  G4double R = fRadiusConst*fCalc->Z13(At);

  static const G4double meanA = 20.;
  if ( At > 20)   // 20.
  {
    R *= (0.8 + 0.2*G4Exp( -((G4double)At - meanA)/meanA) ); 
  }
  else
  {
    R *= (1.0 + 0.1*( 1. - G4Exp( ((G4double)At - meanA)/meanA) ) ); 
  }
  return R;
}

/////////////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::GetNucleusRadiusDE(G4int Z, G4int A)
{
  // algorithm from diffuse-elastic
  static const G4double a11 = 1.26;  // 1.08, 1.16
  static const G4double a12 = 1.19;  // 1.08, 1.16
  static const G4double a13 = 1.12;  // 1.08, 1.16
  static const G4double a2 = 1.1;
  static const G4double a3 = 1.;

  G4double R = CLHEP::fermi;
  // Special rms radii for light nucleii
  if (A < 50)
  {
    if(A == 1)                { return 0.89*R; }// p
    else if(A == 2)           { return 2.13*R; }// d
    else if(Z == 1 && A == 3) { return 1.80*R; }// t
    else if(Z == 2 && A == 3) { return 1.96*R; }// He3
    else if(Z == 2 && A == 4) { return 1.68*R; }// He4
    else if(Z == 3)           { return 2.40*R; }// Li7
    else if(Z == 4)           { return 2.51*R; }// Be9
    else if( 10 < A && A <= 15) { R *= a11*(1. - 1./fCalc->Z23(A)); }
    else if( 15 < A && A <= 20) { R *= a12*(1. - 1./fCalc->Z23(A)); }
    else if( 20 < A && A <= 30) { R *= a13*(1. - 1./fCalc->Z23(A)); }
    else                        { R *= a2; }

    R *= fCalc->Z13(A);
  }
  else
  {
    R *= a3*fCalc->powZ(A, 0.27);
  }
  return R;
}

/////////////////////////////////////////////////////////////////////////////
//
// RMS radii from e-A scattering data

G4double 
G4ComponentGGNuclNuclXsc::GetNucleusRadiusRMS(G4int Z, G4int A)
{
  if     (A == 1)           { return 0.89*fermi; }// p
  else if(A == 2)           { return 2.13*fermi; } // d
  else if(Z == 1 && A == 3) { return 1.80*fermi; }// t

  else if(Z == 2 && A == 3) { return 1.96*fermi; }// He3
  else if(Z == 2 && A == 4) { return 1.68*fermi; }// He4

  else if(Z == 3)           { return 2.40*fermi; }// Li7
  else if(Z == 4)           { return 2.51*fermi; }// Be9

  else                      { return 1.24*fCalc->powZ(A, 0.28 )*fermi; }// A > 9
}

///////////////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::CalculateEcmValue(G4double mp, 
                                                     G4double mt, 
                                                     G4double Plab)
{
  G4double Elab = std::sqrt ( mp * mp + Plab * Plab );
  G4double Ecm  = std::sqrt ( mp * mp + mt * mt + 2 * Elab * mt );
  return Ecm ; // KEcm;
}

///////////////////////////////////////////////////////////////////////////////

G4double G4ComponentGGNuclNuclXsc::CalcMandelstamS(G4double mp, 
                                                   G4double mt, 
                                                   G4double Plab)
{
  G4double Elab = std::sqrt ( mp * mp + Plab * Plab );
  G4double sMand  = mp*mp + mt*mt + 2*Elab*mt ;

  return sMand;
}

///////////////////////////////////////////////////////////////////////////////
