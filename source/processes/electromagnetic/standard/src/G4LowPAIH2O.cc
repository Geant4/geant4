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
//
// 
// G4LowPAIH2O.cc -- class implementation file
//
// GEANT 4 class implementation file
//
// For information related to this code, please, contact
// the Geant4 Collaboration.
//
// R&D: Vladimir.Grichine@cern.ch
//
// History:
//
// 20.07.23 V. Grichine, 1st version 
//

#include "G4LowPAIH2O.hh"
#include "G4LowDataH2O.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4Integrator.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsTable.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4RandomDirection.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4VEmAngularDistribution.hh"


using namespace std;
using namespace CLHEP;

//////////////////////////////////////////////////////////////////
//
// Constructor

G4LowPAIH2O::G4LowPAIH2O( const G4ParticleDefinition* p,
	       const G4String& nam )
    : G4VEmModel(nam) , G4VEmFluctuationModel(nam)
{
  fCof  = fine_structure_const/hbarc/pi;
  fBeta = 0.5;
  fBe2  = fBeta*fBeta;
  fTkin = eV * 1.;
  fBias = 1.; // 
  G4NistManager* man = G4NistManager::Instance();
  fMat = man->FindOrBuildMaterial("G4_WATER");
  // (?)
  fElectronDensity   = fMat->GetElectronDensity();
  fNat = fMat->GetTotNbOfAtomsPerVolume();
  // G4cout<<"fNat = "<<fNat*cm3<<" 1/cm3"<<G4endl;
  fNel = fMat->GetTotNbOfElectPerVolume();
  // G4cout<<"fNel = "<<fNel*cm3<<" 1/cm3"<<G4endl;
  theElectron = G4Electron::Electron();
  theProton   = G4Proton::Proton();
  fParticleChange = nullptr;
  const G4DataVector cuts;
  
  if( IsMaster() )
  {
    InitialiseElementSelectors(p, cuts);
  }
  if(nullptr == fParticleChange)
  {
    fParticleChange = GetParticleChangeForLoss();
  }
  fTotBin = 200; //
  fBinTr  = 100;
  fBmin   = 1.e-3;
  fBmax   = 0.98; // gamma ~ 4 min of dE/dx
  
  fBetaVector = new G4PhysicsLogVector( fBmin, fBmax, fTotBin );

  fWmin = eV  * 0.1; //
  InitRuthELF();
  // BuildPhysicsTable(theProton);
  // BuildPhysicsTable(theElectron);
}

////////////////////////////////////////////////////////////////////////////
//
// Destructor

G4LowPAIH2O::~G4LowPAIH2O()
{
  if(fBetaVector) delete fBetaVector;
  if(fTransferVector) delete fTransferVector;
  if(fPrEnergyTable)
  {
    fPrEnergyTable->clearAndDestroy();
    delete fPrEnergyTable;
  }
  if(fElEnergyTable)
  {
    fElEnergyTable->clearAndDestroy();
    delete fElEnergyTable;
  }
}

//////////////////////////////////////////

void G4LowPAIH2O::Initialise( const G4ParticleDefinition* pd,
			      const G4DataVector&               )
{
  if(nullptr == fParticleChange)
  {
    fParticleChange = GetParticleChangeForLoss();
  }
  // InitRuthELF();
  BuildPhysicsTable( pd ); 
  return;
}

////////////////////////////////

void G4LowPAIH2O::InitialiseLocal(const G4ParticleDefinition* pd,
                                          G4VEmModel* masterModel)
{
  if(nullptr == fParticleChange)
  {
    fParticleChange = GetParticleChangeForLoss();
  }
  SetElementSelectors(masterModel->GetElementSelectors());
  BuildPhysicsTable( pd ); 
  return;
}

////////////////////////////////////////////

void G4LowPAIH2O::Initialize()
{ 
  if( nullptr == fParticleChange )
  {
    fParticleChange = GetParticleChangeForLoss();
  }
  return;
}

///////////////////////////////////////////////////////////
//
// Reading static const arrays to local data vectors for fast searching

void G4LowPAIH2O::InitRuthELF()
{
  G4int i(0), nn(0);
  G4double ee(0.), elf(0.), ruth(0.);

  nn = theBin;

  for( i = 0; i < nn; ++i )
  {
    ee   = theEsum[i];   
    elf  = theELFsum[i];
    ruth = theRuthSum[i];
    fEsum.push_back(ee);
    fELFsum.push_back(elf);
    fRuthSum.push_back(ruth);
  }
  G4cout<<G4endl;
}

//////////////////////////////////////////////////

void G4LowPAIH2O::BuildPhysicsTable(const G4ParticleDefinition* pd )
{
  if( pd == theProton )        BuildPrEnergyTable();
  else if( pd == theElectron ) BuildElEnergyTable();
  else G4cout<<" G4LowPAIH2O::BuildPhysicsTable is not applicable for "
	     <<pd->GetParticleName()<<G4endl;
  return;
}

//////////////////////////////////////////////
//
// Biased dN/dx

G4double G4LowPAIH2O::CrossSectionPerVolume(const G4Material*,
                                 const G4ParticleDefinition* pd,
                                 G4double Tkin,
					    G4double, // cutEnergy,
					    G4double ) // maxEnergy)
{
  G4double dNdx(0.);
  
  if( pd == theProton )        dNdx = GetPrdNdx( Tkin );
  else if( pd == theElectron ) dNdx = GetEldNdx( Tkin );
  else G4cout<<" G4LowPAIH2O::CrossSectionPerVolume is not applicable for "
	     <<pd->GetParticleName()<<G4endl;
  dNdx /= fBias; //

  return dNdx;
}

///////////////////////////

G4double G4LowPAIH2O::CrossSectionPerAtom(  const G4ParticleDefinition* pd,
                                 G4double Tkin,
					    G4double, // Z,
					    G4double, //  A,
                                 G4double cutEnergy,
                                 G4double maxEnergy)
{
  return CrossSectionPerVolume(nullptr, pd, Tkin, cutEnergy, maxEnergy)/fNat;
}

/////////////////////////////

G4double G4LowPAIH2O::ComputeCrossSectionPerElectron(
                                 const G4ParticleDefinition* pd,
                                 G4double Tkin,
                                 G4double cutEnergy,
                                 G4double maxEnergy)
{
  return CrossSectionPerVolume(nullptr, pd, Tkin, cutEnergy, maxEnergy)/fNel;
}
  
////////////////////////////////////////////////

void G4LowPAIH2O::BuildPrEnergyTable()
{
  // G4cout<<"G4LowPAIH2O::BuildPrEnergyTable is called"<<G4endl<<G4endl;
  G4int iBeta(0), iTransfer(0);
  G4double be(0.), be2(0.), gg(1.);
  G4double tp(0.), sum(0.), tmp(0.), mp = proton_mass_c2;

  G4Integrator<G4LowPAIH2O, G4double (G4LowPAIH2O::*)(G4double)> integral;
  fPrEnergyTable = new G4PhysicsTable();

  for( iBeta = 0; iBeta < fTotBin; ++iBeta )  
  {
    be = fBetaVector->GetLowEdgeEnergy(iBeta);
    be2 = be*be;
    SetBe2(be2);
    gg = sqrt( 1./( 1. - be2 ) );
    tp = ( gg - 1.)*mp;
    fWmax = GetProtonTmax(tp);
    if( fWmax <= fWmin ) fWmax = fWmin*1.01;
    fTransferVector = new G4PhysicsLogVector( fWmin, fWmax, fBinTr );
    fPrWmaxVector.push_back(fWmax);
    sum = 0.;
    fTransferVector->PutValue(fBinTr-1,sum);
    
    for( iTransfer = fBinTr-2; iTransfer >= 0; --iTransfer)
    {
      tmp = integral.Legendre10( this, &G4LowPAIH2O::PrPAId2Ndxdw,
				 fTransferVector->GetLowEdgeEnergy(iTransfer),
				 fTransferVector->GetLowEdgeEnergy(iTransfer+1)
			       );
      sum += tmp*fCof/be2;
      // G4cout<<sum*micrometer<<", ";
      fTransferVector->PutValue(iTransfer,sum);
    }
    fPrEnergyTable->insertAt(iBeta, fTransferVector);
  }
  return;
}

////////////////////////////////////////////////

G4double G4LowPAIH2O::GetPrTransfer( G4double Tkin) 
{
  // std::size_t
    G4int iBeta(0), iTransfer(0);
  G4double be(0.), be2(0.), gg(1.), mp = proton_mass_c2;
  G4double rr1(0.), rr2(0.), tr1(0.), tr2(0.), transfer(0.);
  G4double mean1(0.), mean2(0.);
  gg  = 1. + Tkin/mp;
  be2 = 1. - 1./gg/gg;
  if( be2 < 0. ) be2 = 0.;
  be = sqrt(be2);
  
  for( iBeta = 0; iBeta < fTotBin; ++iBeta )  
  {
    if( be <= fBetaVector->GetLowEdgeEnergy(iBeta) ) break;
  }
  if (iBeta == 0 || iBeta == fTotBin-1)
  {
    rr1 = G4UniformRand()*(*(*fPrEnergyTable)(iBeta))(0);

    for( iTransfer = fBinTr-2; iTransfer >= 0; --iTransfer)
    {
      if( (*(*fPrEnergyTable)(iBeta))(iTransfer) >= rr1 ) break;
    }
    if( iTransfer == 0 || iTransfer == fBinTr-2 )
    {
      transfer = (*fPrEnergyTable)(iBeta)->GetLowEdgeEnergy(iTransfer);
    }
    else
    {
      tr1 = (*fPrEnergyTable)(iBeta)->GetLowEdgeEnergy(iTransfer-1); 
      tr2 = (*fPrEnergyTable)(iBeta)->GetLowEdgeEnergy(iTransfer);
      transfer = tr1 + (tr2-tr1)*G4UniformRand();
    }
  }
  else
  {
    rr1 = G4UniformRand()*(*(*fPrEnergyTable)(iBeta-1))(0);

    for( iTransfer = fBinTr-2; iTransfer >= 0; --iTransfer)
    {
      if( (*(*fPrEnergyTable)(iBeta-1))(iTransfer) >= rr1 ) break;
    }
    if( iTransfer == 0 || iTransfer == fBinTr-2 )
    {
      transfer = (*fPrEnergyTable)(iBeta-1)->GetLowEdgeEnergy(iTransfer);
    }
    else
    {
      tr1 = (*fPrEnergyTable)(iBeta-1)->GetLowEdgeEnergy(iTransfer-1); 
      tr2 = (*fPrEnergyTable)(iBeta-1)->GetLowEdgeEnergy(iTransfer);
      mean1 = tr1 + (tr2-tr1)*G4UniformRand();
    }    
    rr2 = G4UniformRand()*(*(*fPrEnergyTable)(iBeta))(0);

    for( iTransfer = fBinTr-2; iTransfer >= 0; --iTransfer)
    {
      if( (*(*fPrEnergyTable)(iBeta))(iTransfer) >= rr2 ) break;
    }
    if( iTransfer == 0 || iTransfer == fBinTr-2 )
    {
      transfer = (*fPrEnergyTable)(iBeta-1)->GetLowEdgeEnergy(iTransfer);
    }
    else
    {
      tr1 = (*fPrEnergyTable)(iBeta)->GetLowEdgeEnergy(iTransfer-1); 
      tr2 = (*fPrEnergyTable)(iBeta)->GetLowEdgeEnergy(iTransfer);
      mean2 = tr1 + (tr2-tr1)*G4UniformRand();
    }
    transfer = mean1 + ( mean2 - mean1 )*G4UniformRand();
  }  
  transfer *= CorrectPrTransfer( Tkin );
  // G4cout<<transfer/eV<<" ";  
  return transfer;
}

/////////////////////////////////////////
//
// Correction for dN/dx -> dE/dx

G4double G4LowPAIH2O::CorrectPrTransfer( G4double Tkin)
{
  // tuning
  G4double Tref = keV * 100.; // 100.; //
  G4double T2   = Tref *  1.; //1.1; //
  G4double T1   = Tref * 1.; // 0.9; // 0.8; //
  G4double rat2 = Tkin/T2;
  G4double rat1 = Tkin/T1;
  G4double yy(0.);
  G4double kk = 0.22; // 0.24; //  
  if( Tkin >= T2 )     yy = kk/pow( rat2, 0.15); //
  else if( Tkin <= T1) yy = kk/pow( rat1, 0.75); //
  else                 yy = kk * 1.; //  1.2; //  1.1; // 

  return yy;
}

/////////////////////////////
//
// Correction for dN/dx -> dE/dx

G4double G4LowPAIH2O::CorrectElTransfer( G4double Tkin)
{
  // tuning
  G4double Tref = eV * 70.; // 100.; //
  G4double T2   = Tref *  1.; //1.1; //
  G4double T1   = Tref * 1.; // 0.9; // 0.8; //
  G4double rat2 = Tkin/T2;
  G4double rat1 = T1/Tkin; // Tkin/T1; //
  G4double yy(0.);
  G4double kk = 0.6; // 0.55; //  0.22; // 
  if( Tkin >= T2 )     yy = kk/pow( rat2, 0.2); // 0.4); //
  else if( Tkin <= T1) yy = kk*pow( rat1, 0.9); // 0.75); //
  else                 yy = kk * 1.; //  1.2; //  1.1; // 

  return yy;
}

////////////////////////////////////////////////

G4double G4LowPAIH2O::GetPrdNdx( G4double Tkin) 
{
  // std::size_t
    G4int iBeta(0);
  G4double be(0.), be2(0.), gg(1.), mp = proton_mass_c2;
  G4double rr1(0.), rr2(0.),  dndx(0.);
  gg  = 1. + Tkin/mp;
  be2 = 1. - 1./gg/gg;
  if(be2 < 0.) be2 = 0.;
  be = sqrt(be2);
  for( iBeta = 0; iBeta < fTotBin; ++iBeta )  
  {
    if( be <= fBetaVector->GetLowEdgeEnergy(iBeta) ) break;
  }
  if (iBeta <= 0 )
  {
    dndx = (*(*fPrEnergyTable)(0))(0);
  }
  else if (iBeta >= fTotBin-1)
  {
    dndx = (*(*fPrEnergyTable)(fTotBin-1))(0);
  }
  else
  {
    rr1 = (*(*fPrEnergyTable)(iBeta-1))(0); 
    rr2 = (*(*fPrEnergyTable)(iBeta))(0);
    dndx = rr1 +(rr2-rr1)*G4UniformRand();
  }
  return dndx;
}

////////////////////////////////////////////////

G4double G4LowPAIH2O::GetPrMFP( G4double Tkin) 
{
  G4double dndx = GetPrdNdx(Tkin);
  G4double mfp(0.);
  if( dndx > 0.) mfp = 1./dndx;
  else           mfp = DBL_MAX;
  return mfp;
}

////////////////////////////////////////////////

G4double G4LowPAIH2O::GetEldNdx( G4double Tkin) 
{
  // std::size_t
    G4int iBeta(0);
  G4double be(0.), be2(0.), gg(1.), me = electron_mass_c2;
  G4double rr1(0.), rr2(0.),  dndx(0.);
  gg  = 1. + Tkin/me;
  be2 = 1. - 1./gg/gg;
  if(be2 < 0.) be2 = 0.;
  be = sqrt(be2);
  for( iBeta = 0; iBeta < fTotBin; ++iBeta )  
  {
    if( be <= fBetaVector->GetLowEdgeEnergy(iBeta) ) break;
  }
  if (iBeta <= 0 )
  {
    dndx = (*(*fElEnergyTable)(0))(0);
  }
  else if (iBeta >= fTotBin-1)
  {
    dndx = (*(*fElEnergyTable)(fTotBin-1))(0);
  }
  else
  {
    rr1  = (*(*fElEnergyTable)(iBeta-1))(0); 
    rr2  = (*(*fElEnergyTable)(iBeta))(0);
    dndx = rr1 + ( rr2 - rr1 )*G4UniformRand();
  }
  return dndx;
}

////////////////////////////////////////////////

G4double G4LowPAIH2O::GetElMFP( G4double Tkin) 
{
  G4double dndx = GetEldNdx(Tkin);
  G4double mfp(0.);
  if( dndx > 0.) mfp = 1./dndx;
  else           mfp = DBL_MAX;  
  return mfp;
}

////////////////////////////////////////////////

void G4LowPAIH2O::BuildElEnergyTable()
{
  // G4cout<<"G4LowPAIH2O::BuildElEnergyTable is called"<<G4endl<<G4endl;
  G4int iBeta(0), iTransfer(0);
  G4double be(0.), be2(0.), gg(1.);
  G4double tp(0.), sum(0.), tmp(0.), me = electron_mass_c2;

  G4Integrator<G4LowPAIH2O, G4double (G4LowPAIH2O::*)(G4double)> integral;
  fElEnergyTable = new G4PhysicsTable();

  for( iBeta = 0; iBeta < fTotBin; ++iBeta )  
  {
    be = fBetaVector->GetLowEdgeEnergy(iBeta);
    be2 = be*be;
    SetBe2(be2);
    gg = sqrt( 1./( 1. - be2 ) );
    tp = ( gg - 1.)*me;
    fWmax = GetElectronTmax(tp);
    if( fWmax <= fWmin ) fWmax = fWmin*1.01;
    fTransferVector = new G4PhysicsLogVector( fWmin, fWmax, fBinTr );

    sum = 0.;
    fTransferVector->PutValue( fBinTr-1, sum );
    
    for( iTransfer = fBinTr-2; iTransfer >= 0; --iTransfer)
    {
      tmp = integral.Legendre10( this, &G4LowPAIH2O::ElPAId2Ndxdw,
				 fTransferVector->GetLowEdgeEnergy(iTransfer),
				 fTransferVector->GetLowEdgeEnergy(iTransfer+1)
			       );
      
      sum += tmp; // *fCof/be2;
      
      fTransferVector->PutValue( iTransfer, sum );
    }
    fElEnergyTable->insertAt( iBeta, fTransferVector );
  }
  return;
}

/////////////////////////////////////////////////////
//
// return d2Ndxdw  (?) collision electron limit

G4double G4LowPAIH2O::ElPAId2Ndxdw( G4double omega )
{
  G4double be2  = GetBe2();
  omega *= 1.7; // tune from p to e-
  
  G4double elf  = GetSumELF(omega);
  G4double ruth = GetSumRuth(omega);
  G4double me   = electron_mass_c2;
  
  G4double d2Ndxdw = elf; // 0.; //

  d2Ndxdw *= log( 2*me*be2/omega );
  
  d2Ndxdw += ruth;
  
  d2Ndxdw *= fCof/be2;

  d2Ndxdw *= 1.3; // norm  to max exp
  
  return d2Ndxdw;
}

////////////////////////////////////////////////

G4double G4LowPAIH2O::GetElTransfer( G4double Tkin) 
{
  // std::size_t
  G4int iBeta, iTransfer;
  G4double be(0.), be2(0.), gg(1.), me = electron_mass_c2;
  G4double rr1(0.), rr2(0.), tr1(0.), tr2(0.), transfer(0.);
  G4double mean1(0.), mean2(0.);
  gg  = 1. + Tkin/me;
  be2 = 1. - 1./gg/gg;
  if( be2 < 0.) be2 = 0.;
  be = sqrt(be2);
  
  for( iBeta = 0; iBeta < fTotBin; ++iBeta )  
  {
    if( be <= fBetaVector->GetLowEdgeEnergy(iBeta) ) break;
  }
  if (iBeta == 0 || iBeta == fTotBin-1)
  {
    rr1 = G4UniformRand()*(*(*fElEnergyTable)(iBeta))(0);

    for( iTransfer = fBinTr-2; iTransfer >= 0; --iTransfer)
    {
      if( (*(*fElEnergyTable)(iBeta))(iTransfer) >= rr1 ) break;
    }
    if( iTransfer == 0 || iTransfer == fBinTr-2 )
    {
      transfer = (*fElEnergyTable)(iBeta)->GetLowEdgeEnergy(iTransfer);
    }
    else
    {
      tr1 = (*fElEnergyTable)(iBeta)->GetLowEdgeEnergy(iTransfer-1); 
      tr2 = (*fElEnergyTable)(iBeta)->GetLowEdgeEnergy(iTransfer);
      transfer = tr1 + (tr2-tr1)*G4UniformRand();
    }
  }
  else
  {
    rr1 = G4UniformRand()*(*(*fElEnergyTable)(iBeta-1))(0);

    for( iTransfer = fBinTr-2; iTransfer >= 0; --iTransfer)
    {
      if( (*(*fElEnergyTable)(iBeta-1))(iTransfer) >= rr1 ) break;
    }
    if( iTransfer == 0 || iTransfer == fBinTr-2 )
    {
      transfer = (*fElEnergyTable)(iBeta-1)->GetLowEdgeEnergy(iTransfer);
    }
    else
    {
      tr1 = (*fElEnergyTable)(iBeta-1)->GetLowEdgeEnergy(iTransfer-1); 
      tr2 = (*fElEnergyTable)(iBeta-1)->GetLowEdgeEnergy(iTransfer);
      mean1 = tr1 + (tr2-tr1)*G4UniformRand();
    }    
    rr2 = G4UniformRand()*(*(*fElEnergyTable)(iBeta))(0);

    for( iTransfer = fBinTr-2; iTransfer >= 0; --iTransfer)
    {
      if( (*(*fElEnergyTable)(iBeta))(iTransfer) >= rr2 ) break;
    }
    if( iTransfer == 0 || iTransfer == fBinTr-2 )
    {
      transfer = (*fElEnergyTable)(iBeta-1)->GetLowEdgeEnergy(iTransfer);
    }
    else
    {
      tr1 = (*fElEnergyTable)(iBeta)->GetLowEdgeEnergy(iTransfer-1); 
      tr2 = (*fElEnergyTable)(iBeta)->GetLowEdgeEnergy(iTransfer);
      mean2 = tr1 + (tr2-tr1)*G4UniformRand();
    }
    transfer = mean1 +(mean2-mean1)*G4UniformRand();
  }
  transfer *= CorrectElTransfer( Tkin );
  if( transfer < 0. ) return 0.;
  // G4cout<<transfer/eV<<" ";    
  return transfer;
}

//////////////////////////////////////////////////
 
G4double G4LowPAIH2O::PrPAId2Ndxdw( G4double omega )
{
  G4double be2  = GetBe2();
  G4double elf  = GetSumELF(omega);
  G4double ruth = GetSumRuth(omega);
  G4double me   = electron_mass_c2;  
  G4double d2Ndxdw = elf;

  d2Ndxdw *= log( 2*me*be2/omega );
  
  d2Ndxdw += ruth;
  
  return d2Ndxdw;
}

///////////////////////////////////////////////

void G4LowPAIH2O::SampleSecondaries( std::vector<G4DynamicParticle*>* vdp,
                                     const G4MaterialCutsCouple*, // matCC,
                                     const G4DynamicParticle* dp,
                                     G4double, // tmin,
                                     G4double) // maxEnergy   )
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  const G4ParticleDefinition* pd = dp->GetDefinition();

  G4double deltaTkin(0.);
  if( pd == theProton )        deltaTkin = GetPrTransfer(kineticEnergy); 
  else if( pd == theElectron ) deltaTkin = GetElTransfer(kineticEnergy);
  else
  {
    G4cout<<" G4LowPAIH2O::SampleSecondaries is not applicable for "
	     <<pd->GetParticleName()<<G4endl;
    return;
  }
  G4ThreeVector direction= dp->GetMomentumDirection();
  fMass = pd->GetPDGMass();
  G4double totalMomentum = sqrt( kineticEnergy*( kineticEnergy + 2.*fMass ) );

  if( !(deltaTkin <= 0.) && !(deltaTkin > 0))
  {
    G4cout<<"G4LowPAIH2O::SampleSecondaries; deltaTkin = "<<deltaTkin/keV
          <<" keV "<< " for Tkin(MeV)= " << kineticEnergy << G4endl;
    return;
  }
  if( deltaTkin < 0.)
  {
    G4cout<<deltaTkin/eV<<"-"<<kineticEnergy/MeV<<", ";
    return;
  }  
  if( kineticEnergy > deltaTkin) kineticEnergy -= deltaTkin;
  else
  {
    deltaTkin = kineticEnergy;
    kineticEnergy = 0.;
  }
  G4double momD = sqrt( deltaTkin*( deltaTkin + 2.*electron_mass_c2 ) );
  G4ThreeVector dirD = G4RandomDirection(); // low e- cloud
  
  G4ThreeVector dir = totalMomentum*direction - dirD*momD;
  direction = dir.unit();
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(direction);

  if( deltaTkin < eV * 1.) // 10.) //
  {
    fParticleChange->ProposeLocalEnergyDeposit(deltaTkin);
  }
  else
  {
    auto deltaRay = new G4DynamicParticle( theElectron, dirD, deltaTkin ); 
    vdp->push_back( deltaRay );
  }
  return;
}

/////////////////////////////////////////

G4double G4LowPAIH2O::SampleFluctuations(const G4MaterialCutsCouple*,
				      const G4DynamicParticle* dp,
					 const G4double, // tcut,
					 const G4double, // tmax,
                                      const G4double length,
					 const G4double) // meanLoss )
{
  G4double eloss(0.), mfp(0.), sumW(0.), sumL(0.), transfer(0.);
  G4double Tkin = dp->GetKineticEnergy();
  const G4ParticleDefinition* pd = dp->GetDefinition();
  
  if( pd != theProton &&  pd != theElectron)   
  {
    eloss = 0.;
    return eloss;
  }
  
  if( pd == theProton )
  {
    do
    {
      mfp      = GetPrMFP(Tkin);
      sumL    += RandExponential::shoot(mfp);
      transfer = GetPrTransfer(Tkin);
      sumW    += transfer;
      // nn++;
      if( Tkin >= transfer ) Tkin    -= transfer;
      else
      {
	transfer = Tkin;
	Tkin     = 0.;
      }
    }
    while( sumL <= length && Tkin >= keV * 0.01 ); //
  }
  else if( pd == theElectron ) // e-
  {
    do
    {
      mfp = GetElMFP(Tkin);
      sumL += RandExponential::shoot(mfp);
      transfer = GetElTransfer(Tkin);
      sumW += transfer;
      Tkin -= transfer;
    }
    while( sumL <= length && Tkin >= 0. );
  }
  eloss = sumW;
  if(eloss < 0.) { eloss = 0.; }
  fParticleChange->SetProposedKineticEnergy( Tkin ); 
  fParticleChange->ProposeLocalEnergyDeposit( eloss );
  // G4cout<<Tkin/MeV<<"/"<<eloss/eV<<"/"<<nn<<", ";
  
  return eloss;
}

/////////////////////////////////////////////

void G4LowPAIH2O::CorrectionsAlongStep(const G4Material*,
				       const G4ParticleDefinition* pd,
				       const G4double Tkin,
				       const G4double,
				       const G4double& length,
				       G4double& eloss)
{
  G4double mfp(0.), sumW(0.);
  
  if     ( pd == theProton )    mfp = GetPrMFP(Tkin);
  else if( pd == theElectron )  mfp = GetElMFP(Tkin);
  else
  {
    eloss = 0.;
    return;
  }
  G4int NN = static_cast<G4int>(RandPoisson::shoot(length/mfp)); 
  
  if( pd == theProton )
  {
    for (G4int nn = 0; nn < NN; ++nn )
    {
      sumW += GetPrTransfer( Tkin );
    }
  }
  else // e-
  {
    for (G4int nn = 0; nn < NN; ++nn )
    {
      sumW += GetElTransfer( Tkin );
    }
  }
  if( sumW < Tkin ) eloss = sumW;
  else              eloss = Tkin;
  // G4cout<<eloss/eV<<".. ";
  fParticleChange->ProposeLocalEnergyDeposit(eloss);
  fParticleChange->SetProposedKineticEnergy(Tkin);
  return;
}

/////////////////////////////////////////////
//
// Transition from secondary e- to radicals
 
G4double G4LowPAIH2O::GetElectronTmax( G4double tt )
{
  fTkin = tt;
  return tt*0.5;

  // algorithm below is doing the same 
  /*  
  G4double tmax(0.), cof(1.);
  G4double mt = eV * 18.; // 32.; // 30.; // 50.; //
  G4double dt = eV * 10.; // 5.; // 20.; // 25.; //

  G4double lim1 = 0.5; // 0.6; //
  G4double lim2 = 1.; // 0.8; //
  G4double dlim = (lim2-lim1)*0.5;

  G4double xx = ( tt - mt )/dt;

  if( xx > 0.) tmax = tt*( lim1 + dlim*exp(-xx) );
  else         tmax = tt*( lim2 - dlim*exp(+xx) );

  dt = mt *0.5; // eV * 38.; // 27.; //
  
  if(tt > mt+dt) tmax = tt*0.5;
  else if (tt < mt-dt) tmax = tt;
  else
  {
    cof = 1. - (tt-mt+dt)*0.25/dt;
    tmax = tt*cof;
  }
  tmax = tt * 0.5; //

  if( tmax > tt ) tmax = tt;
  
  return tmax;
  */
}

/////////////////////////////////////////////////////

G4double G4LowPAIH2O::GetProtonTmax( G4double Tkin )
{
  G4double mp = proton_mass_c2;
  G4double Et = eV * 13.6; // 
  G4double tau      = Tkin/mp;
  G4double gamma    = tau + 1.0;
  G4double bg2      = tau*( tau + 2.0 );
  // G4double beta2    = bg2/(gamma*gamma);
  G4double me = electron_mass_c2;
  G4double rateMass = me/mp;

  if( Tkin < Et ) rateMass = 1.;

  G4double Tmax     = 2.0*me*bg2
                   /( 1. + 2.*gamma*rateMass + rateMass*rateMass );

  return Tmax;
}

//
//
//////////////////////////////////////////////////
