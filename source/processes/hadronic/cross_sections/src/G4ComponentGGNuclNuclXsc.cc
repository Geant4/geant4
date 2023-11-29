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
// 27.05.19 V. Ivantchenko Removed obsolete methods and members 

#include "G4ComponentGGNuclNuclXsc.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadronNucleonXsc.hh"
#include "G4ComponentGGHadronNucleusXsc.hh" 
#include "G4NuclearRadii.hh"
#include "G4Pow.hh"

static const G4double inve = 1./CLHEP::eplus;

G4ComponentGGNuclNuclXsc::G4ComponentGGNuclNuclXsc() 
 : G4VComponentCrossSection("Glauber-Gribov Nucl-nucl"),
   fTotalXsc(0.0), fElasticXsc(0.0), fInelasticXsc(0.0), fProductionXsc(0.0),
   fDiffractionXsc(0.0), fEnergy(0.0), fParticle(nullptr), fZ(0), fA(0)
{
  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theLambda  = G4Lambda::Lambda();
  fHNXsc = new G4HadronNucleonXsc();
  fHadrNucl = new G4ComponentGGHadronNucleusXsc(); 
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
  G4Pow* pG4Pow=G4Pow::GetInstance();
  
  G4int pZ = G4lrint(aParticle->GetPDGCharge()*inve);
  G4int pA = aParticle->GetBaryonNumber();
  G4int pL = aParticle->GetNumberOfLambdasInHypernucleus();
  G4bool pHN = aParticle->IsHypernucleus();
  G4double cHN(0.88);

  // hydrogen
  if(1 == Z && 1 == A) {
    G4double e = kinEnergy*CLHEP::proton_mass_c2/aParticle->GetPDGMass();
    fHadrNucl->ComputeCrossSections( theProton, e, pZ, pA, pL );
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

  G4double tR = G4NuclearRadii::Radius(Z, A);  
  G4double pR = G4NuclearRadii::Radius(pZ, pA);
  
  if(pHN) pR *= std::sqrt( pG4Pow->Z23( pA - pL ) + cHN*pG4Pow->Z23( pL ) )/pG4Pow->Z13(pA);
  
  G4double cB = ComputeCoulombBarier(aParticle, kinEnergy, Z, A, pR, tR);

  if ( cB > 0. ) 
  {
    G4double sigma = (pZ*Z+pN*tN)*fHNXsc->HadronNucleonXscNS(theProton, theProton, pTkin);
    if(pHN) sigma += pL*A*fHNXsc->HadronNucleonXsc(theLambda, theProton, pTkin);
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
