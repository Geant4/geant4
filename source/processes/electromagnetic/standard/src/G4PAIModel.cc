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
// $Id$
//
// -------------------------------------------------------------------
//
// GEANT4 Class
// File name:     G4PAIModel.cc
//
// Author: Vladimir.Grichine@cern.ch on base of Vladimir Ivanchenko model interface
//
// Creation date: 05.10.2003
//
// Modifications:
//
// 17.08.04 V.Grichine, bug fixed for Tkin<=0 in SampleSecondary
// 16.08.04 V.Grichine, bug fixed in massRatio for DEDX, CrossSection, SampleSecondary
// 08.04.05 Major optimisation of internal interfaces (V.Ivantchenko)
// 26.07.09 Fixed logic to work with several materials (V.Ivantchenko)
// 21.11.10 V. Grichine verbose flag for protons and G4PAYySection to check sandia table 
//

#include "G4PAIModel.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4MaterialTable.hh"
#include "G4SandiaTable.hh"
#include "G4OrderedTable.hh"

#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Poisson.hh"
#include "G4Step.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChangeForLoss.hh"

////////////////////////////////////////////////////////////////////////

using namespace std;

G4PAIModel::G4PAIModel(const G4ParticleDefinition* p, const G4String& nam)
  : G4VEmModel(nam),G4VEmFluctuationModel(nam),
  fVerbose(0),
  fLowestGamma(1.005),
  fHighestGamma(10000.),
  fTotBin(200),
  fMeanNumber(20),
  fParticle(0),
  fHighKinEnergy(100.*TeV),
  fTwoln10(2.0*log(10.0)),
  fBg2lim(0.0169),
  fTaulim(8.4146e-3)
{  
  fElectron = G4Electron::Electron();
  fPositron = G4Positron::Positron();

  fPAItransferTable  = 0;
  fPAIdEdxTable      = 0;
  fSandiaPhotoAbsCof = 0;
  fdEdxVector        = 0;
  fLambdaVector      = 0;
  fdNdxCutVector     = 0;
  fParticleEnergyVector = 0;
  fSandiaIntervalNumber = 0;
  fMatIndex = 0;
  fDeltaCutInKinEnergy = 0.0;
  fParticleChange = 0;
  fMaterial = 0;
  fCutCouple = 0;

  if(p) { SetParticle(p); }
  else  { SetParticle(fElectron); }

  isInitialised      = false;
}

////////////////////////////////////////////////////////////////////////////

G4PAIModel::~G4PAIModel()
{
  //  G4cout << "PAI: start destruction" << G4endl;
  if(fParticleEnergyVector) delete fParticleEnergyVector;
  if(fdEdxVector)           delete fdEdxVector ;
  if(fLambdaVector)         delete fLambdaVector;
  if(fdNdxCutVector)        delete fdNdxCutVector;

  if( fPAItransferTable )
    {
      fPAItransferTable->clearAndDestroy();
      delete fPAItransferTable ;
    }
  if( fPAIdEdxTable )
    {
      fPAIdEdxTable->clearAndDestroy();
      delete fPAIdEdxTable ;
    }
  if(fSandiaPhotoAbsCof)
    {
      for(G4int i=0;i<fSandiaIntervalNumber;i++)
	{
	  delete[] fSandiaPhotoAbsCof[i];
	}
      delete[] fSandiaPhotoAbsCof;
    }
  //G4cout << "PAI: end destruction" << G4endl;
}

///////////////////////////////////////////////////////////////////////////////

void G4PAIModel::SetParticle(const G4ParticleDefinition* p)
{
  if(fParticle == p) { return; }
  fParticle = p;
  fMass = fParticle->GetPDGMass();
  fSpin = fParticle->GetPDGSpin();
  G4double q = fParticle->GetPDGCharge()/eplus;
  fChargeSquare = q*q;
  fLowKinEnergy = 0.2*MeV*fMass/proton_mass_c2;
  fRatio = electron_mass_c2/fMass;
  fQc = fMass/fRatio;
  fLowestKineticEnergy  = fMass*(fLowestGamma  - 1.0);
  fHighestKineticEnergy = fMass*(fHighestGamma - 1.0);
}

////////////////////////////////////////////////////////////////////////////

void G4PAIModel::Initialise(const G4ParticleDefinition* p,
			    const G4DataVector&)
{
  if( fVerbose > 0 && p->GetParticleName()=="proton") 
  {
    G4cout<<"G4PAIModel::Initialise for "<<p->GetParticleName()<<G4endl;
    fPAIySection.SetVerbose(1);
  }
  else fPAIySection.SetVerbose(0);

  if(isInitialised) { return; }
  isInitialised = true;

  SetParticle(p);

  fParticleEnergyVector = new G4PhysicsLogVector(fLowestKineticEnergy,
						 fHighestKineticEnergy,
						 fTotBin);

  fParticleChange = GetParticleChangeForLoss();

  // Prepare initialization
  const G4ProductionCutsTable* theCoupleTable =
        G4ProductionCutsTable::GetProductionCutsTable();
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  size_t numOfMat   = G4Material::GetNumberOfMaterials();
  size_t numRegions = fPAIRegionVector.size();

  for(size_t iReg = 0; iReg < numRegions; ++iReg) // region loop
  {
    const G4Region* curReg = fPAIRegionVector[iReg];

    for(size_t jMat = 0; jMat < numOfMat; ++jMat) // region material loop
    {
      fMaterial  = (*theMaterialTable)[jMat];
      fCutCouple = theCoupleTable->GetMaterialCutsCouple( fMaterial, 
					  curReg->GetProductionCuts() );
      //G4cout << "Reg <" <<curReg->GetName() << ">  mat <" 
      //	     << fMaterial->GetName() << ">  fCouple= " 
      //	     << fCutCouple<<"  " << p->GetParticleName() <<G4endl;
      if( fCutCouple ) 
      {
	fMaterialCutsCoupleVector.push_back(fCutCouple);

	fPAItransferTable = new G4PhysicsTable(fTotBin+1);
	fPAIdEdxTable = new G4PhysicsTable(fTotBin+1);

	fDeltaCutInKinEnergy = 
	  (*theCoupleTable->GetEnergyCutsVector(1))[fCutCouple->GetIndex()];
     
	//ComputeSandiaPhotoAbsCof();
	BuildPAIonisationTable();

	fPAIxscBank.push_back(fPAItransferTable);
	fPAIdEdxBank.push_back(fPAIdEdxTable);
	fdEdxTable.push_back(fdEdxVector);

	BuildLambdaVector();
	fdNdxCutTable.push_back(fdNdxCutVector);
	fLambdaTable.push_back(fLambdaVector);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////

void G4PAIModel::InitialiseMe(const G4ParticleDefinition*) 
{}

//////////////////////////////////////////////////////////////////

void G4PAIModel::ComputeSandiaPhotoAbsCof()
{ 
  G4int i, j;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4SandiaTable thisMaterialSandiaTable(fMatIndex) ;
  G4int numberOfElements = 
    (*theMaterialTable)[fMatIndex]->GetNumberOfElements();

  G4int* thisMaterialZ = new G4int[numberOfElements] ;

  for(i=0;i<numberOfElements;i++)  
  {
    thisMaterialZ[i] = 
    (G4int)(*theMaterialTable)[fMatIndex]->GetElement(i)->GetZ() ;
  }  
  fSandiaIntervalNumber = thisMaterialSandiaTable.SandiaIntervals
                           (thisMaterialZ,numberOfElements) ;

  fSandiaIntervalNumber = thisMaterialSandiaTable.SandiaMixing
                           ( thisMaterialZ ,
                             (*theMaterialTable)[fMatIndex]->GetFractionVector() ,
        		     numberOfElements,fSandiaIntervalNumber) ;
   
  fSandiaPhotoAbsCof = new G4double*[fSandiaIntervalNumber] ;

  for(i=0; i<fSandiaIntervalNumber; i++)  
  {
    fSandiaPhotoAbsCof[i] = new G4double[5];
  }
   
  for( i = 0 ; i < fSandiaIntervalNumber ; i++ )
  {
    fSandiaPhotoAbsCof[i][0] = thisMaterialSandiaTable.GetPhotoAbsorpCof(i+1,0); 

    for( j = 1; j < 5 ; j++ )
    {
      fSandiaPhotoAbsCof[i][j] = thisMaterialSandiaTable.GetPhotoAbsorpCof(i+1,j)*
                 (*theMaterialTable)[fMatIndex]->GetDensity() ;
    }
  }
  delete[] thisMaterialZ;
}

////////////////////////////////////////////////////////////////////////////
//
// Build tables for the ionization energy loss
//  the tables are built for MATERIALS
//                           *********

void G4PAIModel::BuildPAIonisationTable()
{
  G4double LowEdgeEnergy , ionloss ;
  G4double tau, Tmax, Tmin, Tkin, deltaLow, /*gamma,*/ bg2 ;

  fdEdxVector = new G4PhysicsLogVector( fLowestKineticEnergy,
					fHighestKineticEnergy,
					fTotBin);
  G4SandiaTable* sandia = fMaterial->GetSandiaTable();

  Tmin = sandia->GetSandiaCofForMaterialPAI(0,0)*keV;
  deltaLow = 100.*eV; // 0.5*eV ;

  for (G4int i = 0 ; i <= fTotBin ; i++)  //The loop for the kinetic energy
  {
    LowEdgeEnergy = fParticleEnergyVector->GetLowEdgeEnergy(i) ;
    tau = LowEdgeEnergy/fMass ;
    //gamma = tau +1. ;
    // G4cout<<"gamma = "<<gamma<<endl ;
    bg2 = tau*( tau + 2. );
    Tmax = MaxSecondaryEnergy(fParticle, LowEdgeEnergy); 
    //    Tmax = std::min(fDeltaCutInKinEnergy, Tmax);
    Tkin = Tmax ;

    if ( Tmax < Tmin + deltaLow )  // low energy safety
      Tkin = Tmin + deltaLow ;

    fPAIySection.Initialize(fMaterial, Tkin, bg2);

    // G4cout<<"ionloss = "<<ionloss*cm/keV<<" keV/cm"<<endl ;
    // G4cout<<"n1 = "<<protonPAI.GetIntegralPAIxSection(1)*cm<<" 1/cm"<<endl ;
    // G4cout<<"protonPAI.GetSplineSize() = "<<
    //    protonPAI.GetSplineSize()<<G4endl<<G4endl ;

    G4int n = fPAIySection.GetSplineSize();
    G4PhysicsFreeVector* transferVector = new G4PhysicsFreeVector(n) ;
    G4PhysicsFreeVector* dEdxVector = new G4PhysicsFreeVector(n);

    for( G4int k = 0 ; k < n; k++ )
    {
      transferVector->PutValue( k ,
                                fPAIySection.GetSplineEnergy(k+1),
                                fPAIySection.GetIntegralPAIySection(k+1) ) ;
      dEdxVector->PutValue( k ,
                                fPAIySection.GetSplineEnergy(k+1),
                                fPAIySection.GetIntegralPAIdEdx(k+1) ) ;
    }
    ionloss = fPAIySection.GetMeanEnergyLoss() ;   //  total <dE/dx>

    if ( ionloss < DBL_MIN) { ionloss = 0.0; }
    fdEdxVector->PutValue(i,ionloss) ;

    fPAItransferTable->insertAt(i,transferVector) ;
    fPAIdEdxTable->insertAt(i,dEdxVector) ;

  }                                        // end of Tkin loop
}

///////////////////////////////////////////////////////////////////////
//
// Build mean free path tables for the delta ray production process
//     tables are built for MATERIALS
//

void G4PAIModel::BuildLambdaVector()
{
  fLambdaVector = new G4PhysicsLogVector( fLowestKineticEnergy,
					  fHighestKineticEnergy,
					  fTotBin                ) ;
  fdNdxCutVector = new G4PhysicsLogVector( fLowestKineticEnergy,
					  fHighestKineticEnergy,
					  fTotBin                ) ;
  if(fVerbose > 1)
  {
    G4cout<<"PAIModel DeltaCutInKineticEnergyNow = "
	  <<fDeltaCutInKinEnergy/keV<<" keV"<<G4endl;
  }
  for (G4int i = 0 ; i <= fTotBin ; i++ )
  {
    G4double dNdxCut = GetdNdxCut(i,fDeltaCutInKinEnergy) ;
    //G4cout << "i= " << i << " x= " << dNdxCut << G4endl;
    if(dNdxCut < 0.0) dNdxCut = 0.0; 
    //    G4double lambda = dNdxCut <= DBL_MIN ? DBL_MAX: 1.0/dNdxCut ;
    G4double lambda = DBL_MAX;
    if(dNdxCut > 0.0) lambda = 1.0/dNdxCut;

    fLambdaVector->PutValue(i, lambda) ;
    fdNdxCutVector->PutValue(i, dNdxCut) ;
  }
}

///////////////////////////////////////////////////////////////////////
//
// Returns integral PAI cross section for energy transfers >= transferCut

G4double  
G4PAIModel::GetdNdxCut( G4int iPlace, G4double transferCut)
{ 
  G4int iTransfer;
  G4double x1, x2, y1, y2, dNdxCut;
  // G4cout<<"iPlace = "<<iPlace<<"; "<<"transferCut = "<<transferCut<<G4endl;
  // G4cout<<"size = "<<G4int((*fPAItransferTable)(iPlace)->GetVectorLength())
  //           <<G4endl;  
  for( iTransfer = 0 ; 
       iTransfer < G4int((*fPAItransferTable)(iPlace)->GetVectorLength()) ; 
       iTransfer++)
  {
    if(transferCut <= (*fPAItransferTable)(iPlace)->GetLowEdgeEnergy(iTransfer))
    {
      break ;
    }
  }  
  if ( iTransfer >= G4int((*fPAItransferTable)(iPlace)->GetVectorLength()) )
  {
      iTransfer = (*fPAItransferTable)(iPlace)->GetVectorLength() - 1 ;
  }
  y1 = (*(*fPAItransferTable)(iPlace))(iTransfer-1) ;
  y2 = (*(*fPAItransferTable)(iPlace))(iTransfer) ;
  if(y1 < 0.0) y1 = 0.0;
  if(y2 < 0.0) y2 = 0.0;
  // G4cout<<"y1 = "<<y1<<"; "<<"y2 = "<<y2<<G4endl;
  x1 = (*fPAItransferTable)(iPlace)->GetLowEdgeEnergy(iTransfer-1) ;
  x2 = (*fPAItransferTable)(iPlace)->GetLowEdgeEnergy(iTransfer) ;
  // G4cout<<"x1 = "<<x1<<"; "<<"x2 = "<<x2<<G4endl;

  if ( y1 == y2 )    dNdxCut = y2 ;
  else
  {
    //  if ( x1 == x2  ) dNdxCut = y1 + (y2 - y1)*G4UniformRand() ;
    //    if ( std::abs(x1-x2) <= eV  ) dNdxCut = y1 + (y2 - y1)*G4UniformRand() ;
    if ( std::abs(x1-x2) <= eV  ) dNdxCut = y1 + (y2 - y1)*0.5 ;
    else             dNdxCut = y1 + (transferCut - x1)*(y2 - y1)/(x2 - x1) ;      
  }
  //  G4cout<<""<<dNdxCut<<G4endl;
  if(dNdxCut < 0.0) dNdxCut = 0.0; 
  return dNdxCut ;
}
///////////////////////////////////////////////////////////////////////
//
// Returns integral dEdx for energy transfers >= transferCut

G4double  
G4PAIModel::GetdEdxCut( G4int iPlace, G4double transferCut)
{ 
  G4int iTransfer;
  G4double x1, x2, y1, y2, dEdxCut;
  //G4cout<<"iPlace = "<<iPlace<<"; "<<"transferCut = "<<transferCut<<G4endl;
  // G4cout<<"size = "<<G4int((*fPAIdEdxTable)(iPlace)->GetVectorLength())
  //           <<G4endl;  
  for( iTransfer = 0 ; 
       iTransfer < G4int((*fPAIdEdxTable)(iPlace)->GetVectorLength()) ; 
       iTransfer++)
  {
    if(transferCut <= (*fPAIdEdxTable)(iPlace)->GetLowEdgeEnergy(iTransfer))
    {
      break ;
    }
  }  
  if ( iTransfer >= G4int((*fPAIdEdxTable)(iPlace)->GetVectorLength()) )
  {
      iTransfer = (*fPAIdEdxTable)(iPlace)->GetVectorLength() - 1 ;
  }
  y1 = (*(*fPAIdEdxTable)(iPlace))(iTransfer-1) ;
  y2 = (*(*fPAIdEdxTable)(iPlace))(iTransfer) ;
  if(y1 < 0.0) y1 = 0.0;
  if(y2 < 0.0) y2 = 0.0;
  //G4cout<<"y1 = "<<y1<<"; "<<"y2 = "<<y2<<G4endl;
  x1 = (*fPAIdEdxTable)(iPlace)->GetLowEdgeEnergy(iTransfer-1) ;
  x2 = (*fPAIdEdxTable)(iPlace)->GetLowEdgeEnergy(iTransfer) ;
  //G4cout<<"x1 = "<<x1<<"; "<<"x2 = "<<x2<<G4endl;

  if ( y1 == y2 )    dEdxCut = y2 ;
  else
  {
    //  if ( x1 == x2  ) dEdxCut = y1 + (y2 - y1)*G4UniformRand() ;
    //    if ( std::abs(x1-x2) <= eV  ) dEdxCut = y1 + (y2 - y1)*G4UniformRand() ;
    if ( std::abs(x1-x2) <= eV  ) dEdxCut = y1 + (y2 - y1)*0.5 ;
    else             dEdxCut = y1 + (transferCut - x1)*(y2 - y1)/(x2 - x1) ;      
  }
  //G4cout<<""<<dEdxCut<<G4endl;
  if(dEdxCut < 0.0) dEdxCut = 0.0; 
  return dEdxCut ;
}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::ComputeDEDXPerVolume(const G4Material*,
					  const G4ParticleDefinition* p,
					  G4double kineticEnergy,
					  G4double cutEnergy)
{
  G4int iTkin,iPlace;
  
  //G4double cut = std::min(MaxSecondaryEnergy(p, kineticEnergy), cutEnergy);
  G4double cut = cutEnergy;

  G4double massRatio  = fMass/p->GetPDGMass();
  G4double scaledTkin = kineticEnergy*massRatio;
  G4double charge     = p->GetPDGCharge();
  G4double charge2    = charge*charge;
  const G4MaterialCutsCouple* matCC = CurrentCouple();

  size_t jMat = 0;
  for(;jMat < fMaterialCutsCoupleVector.size(); ++jMat )
  {
    if( matCC == fMaterialCutsCoupleVector[jMat] ) break;
  }
  //G4cout << "jMat= " << jMat 
  //	 << " jMax= " << fMaterialCutsCoupleVector.size()
  //       << " matCC: " << matCC;
  //if(matCC) G4cout << " mat: " << matCC->GetMaterial()->GetName();
  //  G4cout << G4endl;
  if(jMat == fMaterialCutsCoupleVector.size()) return 0.0;

  fPAIdEdxTable = fPAIdEdxBank[jMat];
  fdEdxVector = fdEdxTable[jMat];
  for(iTkin = 0 ; iTkin <= fTotBin ; iTkin++)
  {
    if(scaledTkin < fParticleEnergyVector->GetLowEdgeEnergy(iTkin)) break ;    
  }
  //G4cout << "E= " << scaledTkin << " iTkin= " << iTkin 
  //	 << "  " << fdEdxVector->GetVectorLength()<<G4endl; 
  iPlace = iTkin - 1;
  if(iPlace < 0) iPlace = 0;
  else if(iPlace >= fTotBin) iPlace = fTotBin-1;
  G4double dEdx = charge2*( (*fdEdxVector)(iPlace) - GetdEdxCut(iPlace,cut) );
  if( dEdx < 0.) dEdx = 0.;
  return dEdx;
}

/////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::CrossSectionPerVolume( const G4Material*,
					    const G4ParticleDefinition* p,
					    G4double kineticEnergy,
					    G4double cutEnergy,
					    G4double maxEnergy  ) 
{
  G4int iTkin,iPlace;
  G4double tmax = std::min(MaxSecondaryEnergy(p, kineticEnergy), maxEnergy);
  if(tmax <= cutEnergy) return 0.0;
  G4double massRatio  = fMass/p->GetPDGMass();
  G4double scaledTkin = kineticEnergy*massRatio;
  G4double charge     = p->GetPDGCharge();
  G4double charge2    = charge*charge, cross, cross1, cross2;
  const G4MaterialCutsCouple* matCC = CurrentCouple();

  size_t jMat = 0;
  for(;jMat < fMaterialCutsCoupleVector.size(); ++jMat )
  {
    if( matCC == fMaterialCutsCoupleVector[jMat] ) break;
  }
  if(jMat == fMaterialCutsCoupleVector.size()) return 0.0;

  fPAItransferTable = fPAIxscBank[jMat];

  for(iTkin = 0 ; iTkin <= fTotBin ; iTkin++)
  {
    if(scaledTkin < fParticleEnergyVector->GetLowEdgeEnergy(iTkin)) break ;    
  }
  iPlace = iTkin - 1;
  if(iPlace < 0) iPlace = 0;
  else if(iPlace >= fTotBin) iPlace = fTotBin-1; 

  //G4cout<<"iPlace = "<<iPlace<<"; tmax = "
  // <<tmax<<"; cutEnergy = "<<cutEnergy<<G4endl;  
  cross1 = GetdNdxCut(iPlace,tmax) ;
  //G4cout<<"cross1 = "<<cross1<<G4endl;  
  cross2 = GetdNdxCut(iPlace,cutEnergy) ;  
  //G4cout<<"cross2 = "<<cross2<<G4endl;  
  cross  = (cross2-cross1)*charge2;
  //G4cout<<"cross = "<<cross<<G4endl;  
  if( cross < DBL_MIN) cross = 0.0;
  //  if( cross2 < DBL_MIN) cross2 = DBL_MIN;

  // return cross2;
  return cross;
}

///////////////////////////////////////////////////////////////////////////
//
// It is analog of PostStepDoIt in terms of secondary electron.
//
 
void G4PAIModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
				   const G4MaterialCutsCouple* matCC,
				   const G4DynamicParticle* dp,
				   G4double tmin,
				   G4double maxEnergy)
{
  size_t jMat = 0;
  for(;jMat < fMaterialCutsCoupleVector.size(); ++jMat )
  {
    if( matCC == fMaterialCutsCoupleVector[jMat] ) break;
  }
  if(jMat == fMaterialCutsCoupleVector.size()) return;

  fPAItransferTable = fPAIxscBank[jMat];
  fdNdxCutVector    = fdNdxCutTable[jMat];

  G4double tmax = std::min(MaxSecondaryKinEnergy(dp), maxEnergy);
  if( tmin >= tmax && fVerbose > 0)
  {
    G4cout<<"G4PAIModel::SampleSecondary: tmin >= tmax "<<G4endl;
  }
  G4ThreeVector direction= dp->GetMomentumDirection();
  G4double particleMass  = dp->GetMass();
  G4double kineticEnergy = dp->GetKineticEnergy();

  G4double massRatio     = fMass/particleMass;
  G4double scaledTkin    = kineticEnergy*massRatio;
  G4double totalEnergy   = kineticEnergy + particleMass;
  G4double pSquare       = kineticEnergy*(totalEnergy+particleMass);
 
  G4double deltaTkin     = GetPostStepTransfer(scaledTkin);

  // G4cout<<"G4PAIModel::SampleSecondaries; deltaKIn = "<<deltaTkin/keV<<" keV "<<G4endl;

  if( deltaTkin <= 0. && fVerbose > 0) 
  {
    G4cout<<"G4PAIModel::SampleSecondary e- deltaTkin = "<<deltaTkin<<G4endl;
  }
  if( deltaTkin <= 0.) return;

  if( deltaTkin > tmax) deltaTkin = tmax;

  G4double deltaTotalMomentum = sqrt(deltaTkin*(deltaTkin + 2. * electron_mass_c2 ));
  G4double totalMomentum      = sqrt(pSquare);
  G4double costheta           = deltaTkin*(totalEnergy + electron_mass_c2)
                                /(deltaTotalMomentum * totalMomentum);

  if( costheta > 0.99999 ) costheta = 0.99999;
  G4double sintheta = 0.0;
  G4double sin2 = 1. - costheta*costheta;
  if( sin2 > 0.) sintheta = sqrt(sin2);

  //  direction of the delta electron
  G4double phi  = twopi*G4UniformRand(); 
  G4double dirx = sintheta*cos(phi), diry = sintheta*sin(phi), dirz = costheta;

  G4ThreeVector deltaDirection(dirx,diry,dirz);
  deltaDirection.rotateUz(direction);
  deltaDirection.unit();

  // primary change
  kineticEnergy -= deltaTkin;
  G4ThreeVector dir = totalMomentum*direction - deltaTotalMomentum*deltaDirection;
  direction = dir.unit();
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(direction);

  // create G4DynamicParticle object for e- delta ray 
  G4DynamicParticle* deltaRay = new G4DynamicParticle;
  deltaRay->SetDefinition(G4Electron::Electron());
  deltaRay->SetKineticEnergy( deltaTkin );  //  !!! trick for last steps /2.0 ???
  deltaRay->SetMomentumDirection(deltaDirection); 

  vdp->push_back(deltaRay);
}


///////////////////////////////////////////////////////////////////////
//
// Returns post step PAI energy transfer > cut electron energy according to passed 
// scaled kinetic energy of particle

G4double  
G4PAIModel::GetPostStepTransfer( G4double scaledTkin )
{  
  // G4cout<<"G4PAIModel::GetPostStepTransfer"<<G4endl ;

  G4int iTkin, iTransfer, iPlace  ;
  G4double transfer = 0.0, position, dNdxCut1, dNdxCut2, E1, E2, W1, W2, W ;

  for(iTkin=0;iTkin<=fTotBin;iTkin++)
  {
    if(scaledTkin < fParticleEnergyVector->GetLowEdgeEnergy(iTkin))  break ;
  }
  iPlace = iTkin - 1 ;
  // G4cout<<"from search, iPlace = "<<iPlace<<G4endl ;
  if(iPlace < 0) iPlace = 0;
  else if(iPlace >= fTotBin) iPlace = fTotBin-1; 
  dNdxCut1 = (*fdNdxCutVector)(iPlace) ;  
  // G4cout<<"dNdxCut1 = "<<dNdxCut1<<G4endl ;


  if(iTkin == fTotBin) // Fermi plato, try from left
  {
      position = dNdxCut1*G4UniformRand() ;

      for( iTransfer = 0;
 iTransfer < G4int((*fPAItransferTable)(iPlace)->GetVectorLength()); iTransfer++ )
      {
        if(position >= (*(*fPAItransferTable)(iPlace))(iTransfer)) break ;
      }
      transfer = GetEnergyTransfer(iPlace,position,iTransfer);
  }
  else
  {
    dNdxCut2 = (*fdNdxCutVector)(iPlace+1) ;  
    // G4cout<<"dNdxCut2 = "<<dNdxCut2<<G4endl ;
    if(iTkin == 0) // Tkin is too small, trying from right only
    {
      position = dNdxCut2*G4UniformRand() ;

      for( iTransfer = 0;
  iTransfer < G4int((*fPAItransferTable)(iPlace+1)->GetVectorLength()); iTransfer++ )
      {
        if(position >= (*(*fPAItransferTable)(iPlace+1))(iTransfer)) break ;
      }
      transfer = GetEnergyTransfer(iPlace+1,position,iTransfer);
    } 
    else // general case: Tkin between two vectors of the material
    {
      E1 = fParticleEnergyVector->GetLowEdgeEnergy(iTkin - 1) ; 
      E2 = fParticleEnergyVector->GetLowEdgeEnergy(iTkin)     ;
      W  = 1.0/(E2 - E1) ;
      W1 = (E2 - scaledTkin)*W ;
      W2 = (scaledTkin - E1)*W ;

      position = ( dNdxCut1*W1 + dNdxCut2*W2 )*G4UniformRand() ;

        // G4cout<<position<<"\t" ;

      G4int iTrMax1, iTrMax2, iTrMax;

      iTrMax1 = G4int((*fPAItransferTable)(iPlace)->GetVectorLength());
      iTrMax2 = G4int((*fPAItransferTable)(iPlace+1)->GetVectorLength());

      if (iTrMax1 >= iTrMax2) iTrMax = iTrMax2;
      else                    iTrMax = iTrMax1;


      for( iTransfer = 0; iTransfer < iTrMax; iTransfer++ )
      {
          if( position >=
          ( (*(*fPAItransferTable)(iPlace))(iTransfer)*W1 +
            (*(*fPAItransferTable)(iPlace+1))(iTransfer)*W2) ) break ;
      }
      transfer = GetEnergyTransfer(iPlace,position,iTransfer);
    }
  } 
  // G4cout<<"PAImodel PostStepTransfer = "<<transfer/keV<<" keV"<<G4endl ; 
  if(transfer < 0.0 ) transfer = 0.0 ;
  // if(transfer < DBL_MIN ) transfer = DBL_MIN;

  return transfer ;
}

///////////////////////////////////////////////////////////////////////
//
// Returns random PAI energy transfer according to passed 
// indexes of particle kinetic

G4double
G4PAIModel::GetEnergyTransfer( G4int iPlace, G4double position, G4int iTransfer )
{ 
  G4int iTransferMax;
  G4double x1, x2, y1, y2, energyTransfer;

  if(iTransfer == 0)
  {
    energyTransfer = (*fPAItransferTable)(iPlace)->GetLowEdgeEnergy(iTransfer);
  }  
  else
  {
    iTransferMax = G4int((*fPAItransferTable)(iPlace)->GetVectorLength());

    if ( iTransfer >= iTransferMax )  iTransfer = iTransferMax - 1;
    
    y1 = (*(*fPAItransferTable)(iPlace))(iTransfer-1);
    y2 = (*(*fPAItransferTable)(iPlace))(iTransfer);

    x1 = (*fPAItransferTable)(iPlace)->GetLowEdgeEnergy(iTransfer-1);
    x2 = (*fPAItransferTable)(iPlace)->GetLowEdgeEnergy(iTransfer);

    if ( x1 == x2 )    energyTransfer = x2;
    else
    {
      if ( y1 == y2  ) energyTransfer = x1 + (x2 - x1)*G4UniformRand();
      else
      {
        energyTransfer = x1 + (position - y1)*(x2 - x1)/(y2 - y1);
      }
    }
  }
  return energyTransfer;
}

///////////////////////////////////////////////////////////////////////

G4double G4PAIModel::SampleFluctuations( const G4Material* material,
                                         const G4DynamicParticle* aParticle,
				               G4double&,
					       G4double& step,
                                               G4double&)
{
  size_t jMat = 0;
  for(;jMat < fMaterialCutsCoupleVector.size(); ++jMat )
  {
    if( material == fMaterialCutsCoupleVector[jMat]->GetMaterial() ) break;
  }
  if(jMat == fMaterialCutsCoupleVector.size()) return 0.0;

  fPAItransferTable = fPAIxscBank[jMat];
  fdNdxCutVector   = fdNdxCutTable[jMat];

  G4int iTkin, iTransfer, iPlace;
  G4long numOfCollisions=0;

  //G4cout<<"G4PAIModel::SampleFluctuations"<<G4endl ;
  //G4cout<<"in: "<<fMaterialCutsCoupleVector[jMat]->GetMaterial()->GetName()<<G4endl;
  G4double loss = 0.0, charge2 ;
  G4double stepSum = 0., stepDelta, lambda, omega; 
  G4double position, E1, E2, W1, W2, W, dNdxCut1, dNdxCut2, meanNumber;
  G4bool numb = true;
  G4double Tkin       = aParticle->GetKineticEnergy() ;
  G4double MassRatio  = fMass/aParticle->GetDefinition()->GetPDGMass() ;
  G4double charge     = aParticle->GetDefinition()->GetPDGCharge() ;
  charge2             = charge*charge ;
  G4double TkinScaled = Tkin*MassRatio ;

  for(iTkin=0;iTkin<=fTotBin;iTkin++)
  {
    if(TkinScaled < fParticleEnergyVector->GetLowEdgeEnergy(iTkin))   break ;
  }
  iPlace = iTkin - 1 ; 
  if(iPlace < 0) iPlace = 0;
  else if(iPlace >= fTotBin) iPlace = fTotBin - 1;
  //G4cout<<"from search, iPlace = "<<iPlace<<G4endl ;
  dNdxCut1 = (*fdNdxCutVector)(iPlace) ;  
  //G4cout<<"dNdxCut1 = "<<dNdxCut1<<G4endl ;

  if(iTkin == fTotBin) // Fermi plato, try from left
  {
    meanNumber =((*(*fPAItransferTable)(iPlace))(0)-dNdxCut1)*step*charge2;
    if(meanNumber < 0.) meanNumber = 0. ;
    //  numOfCollisions = RandPoisson::shoot(meanNumber) ;
    // numOfCollisions = G4Poisson(meanNumber) ;
    if( meanNumber > 0.) lambda = step/meanNumber;
    else                 lambda = DBL_MAX;
    while(numb)
    {
     stepDelta = CLHEP::RandExponential::shoot(lambda);
     stepSum += stepDelta;
     if(stepSum >= step) break;
     numOfCollisions++;
    }   
    //G4cout<<"##1 numOfCollisions = "<<numOfCollisions<<G4endl ;

    while(numOfCollisions)
    {
      position = dNdxCut1+
                 ((*(*fPAItransferTable)(iPlace))(0)-dNdxCut1)*G4UniformRand() ;

      for( iTransfer = 0;
   iTransfer < G4int((*fPAItransferTable)(iPlace)->GetVectorLength()); iTransfer++ )
      {
        if(position >= (*(*fPAItransferTable)(iPlace))(iTransfer)) break ;
      }
      omega = GetEnergyTransfer(iPlace,position,iTransfer);
      //G4cout<<"G4PAIModel::SampleFluctuations, omega = "<<omega/keV<<" keV; "<<"\t";
      loss += omega;
      numOfCollisions-- ;
    }
  }
  else
  {
    dNdxCut2 = (*fdNdxCutVector)(iPlace+1) ; 
    //G4cout<<"dNdxCut2 = "<<dNdxCut2<< " iTkin= "<<iTkin<<" iPlace= "<<iPlace<<G4endl;
 
    if(iTkin == 0) // Tkin is too small, trying from right only
    {
      meanNumber =((*(*fPAItransferTable)(iPlace+1))(0)-dNdxCut2)*step*charge2;
      if( meanNumber < 0. ) meanNumber = 0. ;
      //  numOfCollisions = CLHEP::RandPoisson::shoot(meanNumber) ;
      //  numOfCollisions = G4Poisson(meanNumber) ;
      if( meanNumber > 0.) lambda = step/meanNumber;
      else                 lambda = DBL_MAX;
      while(numb)
	{
	  stepDelta = CLHEP::RandExponential::shoot(lambda);
	  stepSum += stepDelta;
	  if(stepSum >= step) break;
	  numOfCollisions++;
	}   

      //G4cout<<"##2 numOfCollisions = "<<numOfCollisions<<G4endl ;

      while(numOfCollisions)
      {
        position = dNdxCut2+
                   ((*(*fPAItransferTable)(iPlace+1))(0)-dNdxCut2)*G4UniformRand();
   
        for( iTransfer = 0;
   iTransfer < G4int((*fPAItransferTable)(iPlace+1)->GetVectorLength()); iTransfer++ )
        {
          if(position >= (*(*fPAItransferTable)(iPlace+1))(iTransfer)) break ;
        }
        omega = GetEnergyTransfer(iPlace,position,iTransfer);
        //G4cout<<omega/keV<<"\t";
        loss += omega;
        numOfCollisions-- ;
      }
    } 
    else // general case: Tkin between two vectors of the material
    {
      E1 = fParticleEnergyVector->GetLowEdgeEnergy(iTkin - 1) ; 
      E2 = fParticleEnergyVector->GetLowEdgeEnergy(iTkin)     ;
       W = 1.0/(E2 - E1) ;
      W1 = (E2 - TkinScaled)*W ;
      W2 = (TkinScaled - E1)*W ;

      //G4cout << fPAItransferTable->size() << G4endl;
      //G4cout<<"(*(*fPAItransferTable)(iPlace))(0) = "<<
      //   (*(*fPAItransferTable)(iPlace))(0)<<G4endl ;
      //G4cout<<"(*(*fPAItransferTable)(iPlace+1))(0) = "<<
      //     (*(*fPAItransferTable)(iPlace+1))(0)<<G4endl ;

      meanNumber=( ((*(*fPAItransferTable)(iPlace))(0)-dNdxCut1)*W1 + 
		   ((*(*fPAItransferTable)(iPlace+1))(0)-dNdxCut2)*W2 )*step*charge2;
      if(meanNumber<0.0) meanNumber = 0.0;
      //  numOfCollisions = RandPoisson::shoot(meanNumber) ;
      // numOfCollisions = G4Poisson(meanNumber) ;
      if( meanNumber > 0.) lambda = step/meanNumber;
      else                 lambda = DBL_MAX;
      while(numb)
	{
	  stepDelta = CLHEP::RandExponential::shoot(lambda);
	  stepSum += stepDelta;
	  if(stepSum >= step) break;
	  numOfCollisions++;
	}   

      //G4cout<<"##3 numOfCollisions = "<<numOfCollisions<<endl ;

      while(numOfCollisions)
      {
        position = dNdxCut1*W1 + dNdxCut2*W2 +
                 ( ( (*(*fPAItransferTable)(iPlace))(0)-dNdxCut1 )*W1 + 
                    dNdxCut2+
                  ( (*(*fPAItransferTable)(iPlace+1))(0)-dNdxCut2 )*W2 )*G4UniformRand();

        // G4cout<<position<<"\t" ;

        for( iTransfer = 0;
    iTransfer < G4int((*fPAItransferTable)(iPlace)->GetVectorLength()); iTransfer++ )
        {
          if( position >=
          ( (*(*fPAItransferTable)(iPlace))(iTransfer)*W1 + 
            (*(*fPAItransferTable)(iPlace+1))(iTransfer)*W2) )
          {
	      break ;
	  }
        }
	omega =  GetEnergyTransfer(iPlace,position,iTransfer);
	//  G4cout<<omega/keV<<"\t";
        loss += omega;
        numOfCollisions-- ;    
      }
    }
  } 
  //G4cout<<"PAIModel AlongStepLoss = "<<loss/keV<<" keV, on step = "
  //  <<step/mm<<" mm"<<G4endl ; 
  if(loss > Tkin) loss=Tkin;
  if(loss < 0.  ) loss = 0.;
  return loss ;

}

//////////////////////////////////////////////////////////////////////
//
// Returns the statistical estimation of the energy loss distribution variance
//


G4double G4PAIModel::Dispersion( const G4Material* material, 
                                 const G4DynamicParticle* aParticle,
 				       G4double& tmax, 
			               G4double& step       )
{
  G4double loss, sumLoss=0., sumLoss2=0., sigma2, meanLoss=0.;
  for(G4int i = 0 ; i < fMeanNumber; i++)
  {
    loss      = SampleFluctuations(material,aParticle,tmax,step,meanLoss);
    sumLoss  += loss;
    sumLoss2 += loss*loss;
  }
  meanLoss = sumLoss/fMeanNumber;
  sigma2   = meanLoss*meanLoss + (sumLoss2-2*sumLoss*meanLoss)/fMeanNumber;
  return sigma2;
}

/////////////////////////////////////////////////////////////////////

G4double G4PAIModel::MaxSecondaryEnergy( const G4ParticleDefinition* p,
					 G4double kinEnergy) 
{
  G4double tmax = kinEnergy;
  if(p == fElectron) tmax *= 0.5;
  else if(p != fPositron) { 
    G4double mass = p->GetPDGMass();
    G4double ratio= electron_mass_c2/mass;
    G4double gamma= kinEnergy/mass + 1.0;
    tmax = 2.0*electron_mass_c2*(gamma*gamma - 1.) /
                  (1. + 2.0*gamma*ratio + ratio*ratio);
  }
  return tmax;
}

///////////////////////////////////////////////////////////////

void G4PAIModel::DefineForRegion(const G4Region* r) 
{
  fPAIRegionVector.push_back(r);
}

//
//
/////////////////////////////////////////////////






