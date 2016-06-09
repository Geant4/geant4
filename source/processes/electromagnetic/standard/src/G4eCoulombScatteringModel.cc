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
// $Id: G4eCoulombScatteringModel.cc,v 1.11 2007/05/22 17:34:36 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eCoulombScatteringModel
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 22.08.2005
//
// Modifications:
// 01.08.06 V.Ivanchenko extend upper limit of table to TeV and review the
//          logic of building - only elements from G4ElementTable
// 08.08.06 V.Ivanchenko build internal table in ekin scale, introduce faclim
// 19.08.06 V.Ivanchenko add inline function ScreeningParameter 
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eCoulombScatteringModel.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForGamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eCoulombScatteringModel::G4eCoulombScatteringModel(
  G4double thetaMin, G4double thetaMax, G4bool build, 
  G4double tlim, const G4String& nam)
  : G4VEmModel(nam),
    cosThetaMin(cos(thetaMin)),
    cosThetaMax(cos(thetaMax)),
    q2Limit(tlim),
    theCrossSectionTable(0),
    lowKEnergy(keV),
    highKEnergy(TeV),
    alpha2(fine_structure_const*fine_structure_const),
    faclim(100.0),
    nbins(12),
    nmax(100),
    buildTable(build),
    isInitialised(false)
{
  a0 = alpha2*electron_mass_c2*electron_mass_c2/(0.885*0.885);
  G4double p0 = electron_mass_c2*classic_electr_radius;
  coeff = twopi*p0*p0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eCoulombScatteringModel::~G4eCoulombScatteringModel()
{
  if(theCrossSectionTable) {
    theCrossSectionTable->clearAndDestroy();
    delete theCrossSectionTable;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eCoulombScatteringModel::Initialise(const G4ParticleDefinition* p,
					   const G4DataVector&)
{
  if(isInitialised) return;
  isInitialised = true;

  if(pParticleChange)
    fParticleChange = 
      reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
  else
    fParticleChange = new G4ParticleChangeForGamma();

  if(!buildTable || p->GetParticleName() == "GenericIon") return;

  // Compute cross section multiplied by Ptot^2*beta^2
  theCrossSectionTable = new G4PhysicsTable();
  G4PhysicsLogVector* ptrVector;
  G4double e, value;
  nbins = 2*G4int(log10(highKEnergy/lowKEnergy));

  const  G4ElementTable* elmt = G4Element::GetElementTable();
  size_t nelm =  G4Element::GetNumberOfElements();

  for(size_t j=0; j<nelm; j++) { 

    ptrVector  = new G4PhysicsLogVector(lowKEnergy, highKEnergy, nbins);
    const G4Element* elm = (*elmt)[j]; 
    G4double Z =  elm->GetZ();
    index[G4int(Z)] = j;
    for(G4int i=0; i<=nbins; i++) {
      e     = ptrVector->GetLowEdgeEnergy( i ) ;
      value = CalculateCrossSectionPerAtom(p, e, Z);  
      ptrVector->PutValue( i, log(value) );
    }

    theCrossSectionTable->insert(ptrVector);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eCoulombScatteringModel::CalculateCrossSectionPerAtom(
		             const G4ParticleDefinition* p,      
			     G4double kinEnergy, 
			     G4double Z)
{
  G4double cross = 0.0;
  G4double m     = p->GetPDGMass();
  G4double tkin  = std::max(keV, kinEnergy);
  G4double mom2  = tkin*(tkin + 2.0*m);
  G4double costm = std::max(cosThetaMax, 1.0 - 0.5*q2Limit/mom2);
  if(costm < cosThetaMin) {
    G4double q        = p->GetPDGCharge()/eplus;
    G4double q2       = q*q;
    G4double invbeta2 = 1.0 +  m*m/mom2;
    G4double A = ScreeningParameter(Z, q2, mom2, invbeta2);
    G4double a = 2.0*A + 1.0;
    cross = coeff*Z*(Z + 1.0)*q2*invbeta2*(cosThetaMin - costm)/
      ((a - cosThetaMin)*(a - costm)*mom2);
    /*
    if(Z == 13 || Z == 79) {
      G4cout << "## e= " << kinEnergy << "  beta= " << sqrt (1.0/invbeta2)
	     <<"  Z= " << Z 
	     << " sig(bn)= " << cross/barn 
	     << "  cosMax= " <<  costm 
	     << "  cosMin= " <<  cosThetaMin 
	     << G4endl;
      G4double atommass = 27.0;
      if(Z == 79) atommass = 197.0;
      G4double u0 = 1.e+6*atommass*cm2/(cross*Avogadro);
      G4double u1 = 0.5*u0/( A* ( (1.0 + A)*log(1.0 + 1.0/A) -1.0 ) );
      G4cout << "  l0= " << u0 << "  l1= " << u1 
	     << "   A= " << A << G4endl;
    }
    */
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eCoulombScatteringModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
						  const G4MaterialCutsCouple* couple,
						  const G4DynamicParticle* dp,
						  G4double,
						  G4double)
{
  const G4Material* aMaterial = couple->GetMaterial();
  const G4ParticleDefinition* p = dp->GetDefinition();
  
  G4double mass      = dp->GetMass();
  G4double kinEnergy = dp->GetKineticEnergy();
  G4double mom2      = kinEnergy*(kinEnergy + 2.0*mass);

  const G4Element* elm = SelectRandomAtom(aMaterial, p, kinEnergy);
  G4double Z  = elm->GetZ();
  G4double q  = p->GetPDGCharge()/eplus;
  G4double q2 = q*q;

  G4double invbeta2  = 1.0 + mass*mass/mom2;
  G4double a = 2.*ScreeningParameter(Z, q2, mom2, invbeta2);

  G4double costm = std::max(cosThetaMax, 1.0 - 0.5*q2Limit/mom2);
  if(costm >= cosThetaMin) return; 

  G4double x   = G4UniformRand();
  G4double y   = (a + 1.0 - cosThetaMin)/(cosThetaMin - costm);
  G4double st2 = 0.5*(y*(1.0 - costm) - a*x)/(y + x); 
  if(st2 < 0.0 || st2 > 1.0) {
    G4cout << "G4eCoulombScatteringModel::SampleSecondaries WARNING st2= " 
	   << st2 << G4endl;
    st2 = 0.0;
  }

  G4double tet = 2.0*asin(sqrt(st2));
  G4double cost= cos(tet);
  G4double sint= sin(tet);

  G4double phi  = twopi * G4UniformRand();

  G4ThreeVector direction = dp->GetMomentumDirection(); 
  G4ThreeVector newDirection(cos(phi)*sint,sin(phi)*sint,cost);
  newDirection.rotateUz(direction);   

  fParticleChange->ProposeMomentumDirection(newDirection);   
 
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


