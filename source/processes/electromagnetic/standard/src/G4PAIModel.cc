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
// File name:     G4PAIModel.cc
//
// Author: Vladimir.Grichine@cern.ch on base of Vladimir Ivanchenko code
//
// Creation date: 05.10.2003
//
// Modifications:
//
// 17.08.04 V.Grichine, bug fixed for Tkin<=0 in SampleSecondary
// 16.08.04 V.Grichine, bug fixed in massRatio for DEDX, CrossSection, SampleSecondary
// 08.04.05 Major optimisation of internal interfaces (V.Ivantchenko)
//

#include "G4Region.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4MaterialTable.hh"
#include "G4SandiaTable.hh"
#include "G4PAIxSection.hh"

#include "G4PAIModel.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
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
  fLowKinEnergy(2.0*MeV),
  fTwoln10(2.0*log(10.0)),
  fBg2lim(0.0169),
  fTaulim(8.4146e-3)
{
  if(p) SetParticle(p);
  
  // fLowestKineticEnergy  = LowEnergyLimit();
  // fHighestKineticEnergy = HighEnergyLimit();

  fLowestKineticEnergy  = fMass*fLowestGamma;
  fHighestKineticEnergy = fMass*fHighestGamma;

  fParticleEnergyVector = new G4PhysicsLogVector(fLowestKineticEnergy,
						 fHighestKineticEnergy,
						 fTotBin                );
  fPAItransferTable  = 0;
  fPAIdEdxTable      = 0;
  fSandiaPhotoAbsCof = 0;
  fdEdxVector        = 0;
  fLambdaVector      = 0;
  fdNdxCutVector     = 0;
}

////////////////////////////////////////////////////////////////////////////

G4PAIModel::~G4PAIModel()
{
  if(fParticleEnergyVector) delete fParticleEnergyVector;
  if(fdEdxVector)         delete fdEdxVector ;
  if ( fLambdaVector)     delete fLambdaVector;
  if ( fdNdxCutVector)    delete fdNdxCutVector;

  if( fPAItransferTable )
  {
        fPAItransferTable->clearAndDestroy();
        delete fPAItransferTable ;
  }
  if(fSandiaPhotoAbsCof)
  {
    for(G4int i=0;i<fSandiaIntervalNumber;i++)
    {
        delete[] fSandiaPhotoAbsCof[i];
    }
    delete[] fSandiaPhotoAbsCof;
  }
}

///////////////////////////////////////////////////////////////////////////////

void G4PAIModel::SetParticle(const G4ParticleDefinition* p)
{
  fParticle = p;
  fMass = fParticle->GetPDGMass();
  fSpin = fParticle->GetPDGSpin();
  G4double q = fParticle->GetPDGCharge()/eplus;
  fChargeSquare = q*q;
  fLowKinEnergy *= fMass/proton_mass_c2;
  fRatio = electron_mass_c2/fMass;
  fQc = fMass/fRatio;
}

////////////////////////////////////////////////////////////////////////////

void G4PAIModel::Initialise(const G4ParticleDefinition* p,
			    const G4DataVector&)
{
  if(!fParticle) SetParticle(p);

  if(pParticleChange)
    fParticleChange = reinterpret_cast<G4ParticleChangeForLoss*>(pParticleChange);
  else
    fParticleChange = new G4ParticleChangeForLoss();

  const G4ProductionCutsTable* theCoupleTable =
        G4ProductionCutsTable::GetProductionCutsTable();

  for(size_t iReg = 0; iReg < fPAIRegionVector.size();++iReg) // region loop
  {
    const G4Region* curReg = fPAIRegionVector[iReg];

    // (*fPAIRegionVector[iRegion])

    vector<G4Material*>::const_iterator matIter = curReg->GetMaterialIterator();
    size_t jMat; 
    size_t numOfMat = curReg->GetNumberOfMaterials();

    //  for(size_t jMat = 0; jMat < curReg->GetNumberOfMaterials();++jMat){}
    const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    size_t numberOfMat = G4Material::GetNumberOfMaterials();

    for(jMat = 0 ; jMat < numOfMat; ++jMat) // region material loop
    {
      const G4MaterialCutsCouple* matCouple = theCoupleTable->
      GetMaterialCutsCouple( *matIter, curReg->GetProductionCuts() );
      fMaterialCutsCoupleVector.push_back(matCouple);

      size_t iMatGlob;
      for(iMatGlob = 0 ; iMatGlob < numberOfMat ; iMatGlob++ )
      {
        if( *matIter == (*theMaterialTable)[iMatGlob]) break ;
      }
      fMatIndex = iMatGlob;

      ComputeSandiaPhotoAbsCof();
      BuildPAIonisationTable(p);

      fPAIxscBank.push_back(fPAItransferTable);
      fPAIdEdxBank.push_back(fPAIdEdxTable);
      fdEdxTable.push_back(fdEdxVector);

      BuildLambdaVector(matCouple);
      fdNdxCutTable.push_back(fdNdxCutVector);
      fLambdaTable.push_back(fLambdaVector);


      matIter++;
    }
  }
}

//////////////////////////////////////////////////////////////////

void G4PAIModel::ComputeSandiaPhotoAbsCof()
{
  G4int i, j, numberOfElements ;
  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4SandiaTable thisMaterialSandiaTable(fMatIndex) ;
  numberOfElements = (*theMaterialTable)[fMatIndex]->
                                              GetNumberOfElements();
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

  for(i=0;i<fSandiaIntervalNumber;i++)  fSandiaPhotoAbsCof[i] = new G4double[5] ;
   
  for( i = 0 ; i < fSandiaIntervalNumber ; i++ )
  {
    fSandiaPhotoAbsCof[i][0] = thisMaterialSandiaTable.GetPhotoAbsorpCof(i+1,0) ; 

    for( j = 1; j < 5 ; j++ )
    {
      fSandiaPhotoAbsCof[i][j] = thisMaterialSandiaTable.
	                              GetPhotoAbsorpCof(i+1,j)*
                 (*theMaterialTable)[fMatIndex]->GetDensity() ;
    }
  }
  // delete[] thisMaterialZ ;
}

////////////////////////////////////////////////////////////////////////////
//
// Build tables for the ionization energy loss
//  the tables are built for MATERIALS
//                           *********

void
G4PAIModel::BuildPAIonisationTable(const G4ParticleDefinition* p)
{
  G4double LowEdgeEnergy , ionloss ;
  G4double massRatio, tau, Tmax, Tmin, Tkin, deltaLow, gamma, bg2 ;
  /*
  if( fPAItransferTable )
  {
     fPAItransferTable->clearAndDestroy() ;
     delete fPAItransferTable ;
  }
  */
  fPAItransferTable = new G4PhysicsTable(fTotBin);
  /*
  if( fPAIdEdxTable )
  {
     fPAIdEdxTable->clearAndDestroy() ;
     delete fPAIdEdxTable ;
  }
  */
  fPAIdEdxTable = new G4PhysicsTable(fTotBin);

  //  if(fdEdxVector) delete fdEdxVector ;
  fdEdxVector = new G4PhysicsLogVector( fLowestKineticEnergy,
					 fHighestKineticEnergy,
					 fTotBin               ) ;
  Tmin     = fSandiaPhotoAbsCof[0][0] ;      // low energy Sandia interval
  deltaLow = 100.*eV; // 0.5*eV ;

  for (G4int i = 0 ; i < fTotBin ; i++)  //The loop for the kinetic energy
  {
    LowEdgeEnergy = fParticleEnergyVector->GetLowEdgeEnergy(i) ;
    tau = LowEdgeEnergy/fMass ;
    //    if(tau < 0.01)  tau = 0.01 ;
    gamma = tau +1. ;
    // G4cout<<"gamma = "<<gamma<<endl ;
    bg2 = tau*( tau + 2. );

    if(p->GetParticleName() == "e-")
    {
      Tmax = 0.5*LowEdgeEnergy;
    }
    else if(p->GetParticleName() == "e+")
    {
      Tmax = LowEdgeEnergy;   // Unclear??
    }
    else
    {
      massRatio = electron_mass_c2/fMass ;
      Tmax = 2.*electron_mass_c2*bg2/(1. + 2.*gamma*massRatio+massRatio*massRatio);
    }


    // G4cout<<"proton Tkin = "<<LowEdgeEnergy/MeV<<" MeV"
    // <<" Tmax = "<<Tmax/MeV<<" MeV"<<G4endl;
    // Tkin = DeltaCutInKineticEnergyNow ;

    // if ( DeltaCutInKineticEnergyNow > Tmax)         // was <
    {
      Tkin = Tmax ;
    }
    if ( Tkin < Tmin + deltaLow )  // low energy safety
    {
      Tkin = Tmin + deltaLow ;
    }
    G4PAIxSection protonPAI( fMatIndex,
                             Tkin,
                             bg2,
                             fSandiaPhotoAbsCof,
                             fSandiaIntervalNumber  ) ;


    // G4cout<<"ionloss = "<<ionloss*cm/keV<<" keV/cm"<<endl ;
    // G4cout<<"n1 = "<<protonPAI.GetIntegralPAIxSection(1)*cm<<" 1/cm"<<endl ;
    // G4cout<<"protonPAI.GetSplineSize() = "<<
    //    protonPAI.GetSplineSize()<<G4endl<<G4endl ;

    G4PhysicsFreeVector* transferVector = new
                             G4PhysicsFreeVector(protonPAI.GetSplineSize()) ;
    G4PhysicsFreeVector* dEdxVector = new
                             G4PhysicsFreeVector(protonPAI.GetSplineSize()) ;

    for( G4int k = 0 ; k < protonPAI.GetSplineSize() ; k++ )
    {
      transferVector->PutValue( k ,
                                protonPAI.GetSplineEnergy(k+1),
                                protonPAI.GetIntegralPAIxSection(k+1) ) ;
      dEdxVector->PutValue( k ,
                                protonPAI.GetSplineEnergy(k+1),
                                protonPAI.GetIntegralPAIdEdx(k+1) ) ;
    }
    ionloss = protonPAI.GetMeanEnergyLoss() ;   //  total <dE/dx>
    if ( ionloss <= 0.)  ionloss = DBL_MIN ;
    fdEdxVector->PutValue(i,ionloss) ;

    fPAItransferTable->insertAt(i,transferVector) ;
    fPAIdEdxTable->insertAt(i,dEdxVector) ;

    // delete[] transferVector ;
    }                                        // end of Tkin loop
   //  theLossTable->insert(fdEdxVector);
                                              // end of material loop
   // G4cout<<"G4PAIonisation::BuildPAIonisationTable() have been called"<<G4endl ;
   // G4cout<<"G4PAIonisation::BuildLossTable() have been called"<<G4endl ;
}

///////////////////////////////////////////////////////////////////////
//
// Build mean free path tables for the delta ray production process
//     tables are built for MATERIALS
//

void
G4PAIModel::BuildLambdaVector(const G4MaterialCutsCouple* matCutsCouple)
{
  G4int i ;
  G4double dNdxCut, lambda;

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();

  size_t numOfCouples = theCoupleTable->GetTableSize();
  size_t jMatCC;

  for (jMatCC = 0 ; jMatCC < numOfCouples ; jMatCC++ )
  {
    if( matCutsCouple == theCoupleTable->GetMaterialCutsCouple(jMatCC) ) break;
  }
  if( jMatCC == numOfCouples && jMatCC > 0 ) jMatCC--;

  const vector<G4double>*  deltaCutInKineticEnergy = theCoupleTable->
                                GetEnergyCutsVector(idxG4ElectronCut);

  if (fLambdaVector)   delete fLambdaVector;
  if (fdNdxCutVector)  delete fdNdxCutVector;

  fLambdaVector = new G4PhysicsLogVector( fLowestKineticEnergy,
					  fHighestKineticEnergy,
					  fTotBin                ) ;
  fdNdxCutVector = new G4PhysicsLogVector( fLowestKineticEnergy,
					  fHighestKineticEnergy,
					  fTotBin                ) ;
  G4double deltaCutInKineticEnergyNow = (*deltaCutInKineticEnergy)[jMatCC] ;
  if(fVerbose)
  {
  G4cout<<"PAIModel DeltaCutInKineticEnergyNow = "
        <<deltaCutInKineticEnergyNow/keV<<" keV"<<G4endl;
  }
  for ( i = 0 ; i < fTotBin ; i++ )
  {
    dNdxCut = GetdNdxCut(i,deltaCutInKineticEnergyNow) ;
    lambda = dNdxCut <= DBL_MIN ? DBL_MAX: 1.0/dNdxCut ;

    if (lambda <= 1000*kCarTolerance) lambda = 1000*kCarTolerance ; // Mmm ???

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
  // G4cout<<"y1 = "<<y1<<"; "<<"y2 = "<<y2<<G4endl;
  x1 = (*fPAItransferTable)(iPlace)->GetLowEdgeEnergy(iTransfer-1) ;
  x2 = (*fPAItransferTable)(iPlace)->GetLowEdgeEnergy(iTransfer) ;
  // G4cout<<"x1 = "<<x1<<"; "<<"x2 = "<<x2<<G4endl;

  if ( y1 == y2 )    dNdxCut = y2 ;
  else
  {
    //  if ( x1 == x2  ) dNdxCut = y1 + (y2 - y1)*G4UniformRand() ;
    if ( std::abs(x1-x2) <= eV  ) dNdxCut = y1 + (y2 - y1)*G4UniformRand() ;
    else             dNdxCut = y1 + (transferCut - x1)*(y2 - y1)/(x2 - x1) ;      
  }
  //  G4cout<<""<<dNdxCut<<G4endl;
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
  // G4cout<<"iPlace = "<<iPlace<<"; "<<"transferCut = "<<transferCut<<G4endl;
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
  // G4cout<<"y1 = "<<y1<<"; "<<"y2 = "<<y2<<G4endl;
  x1 = (*fPAIdEdxTable)(iPlace)->GetLowEdgeEnergy(iTransfer-1) ;
  x2 = (*fPAIdEdxTable)(iPlace)->GetLowEdgeEnergy(iTransfer) ;
  // G4cout<<"x1 = "<<x1<<"; "<<"x2 = "<<x2<<G4endl;

  if ( y1 == y2 )    dEdxCut = y2 ;
  else
  {
    //  if ( x1 == x2  ) dEdxCut = y1 + (y2 - y1)*G4UniformRand() ;
    if ( std::abs(x1-x2) <= eV  ) dEdxCut = y1 + (y2 - y1)*G4UniformRand() ;
    else             dEdxCut = y1 + (transferCut - x1)*(y2 - y1)/(x2 - x1) ;      
  }
  //  G4cout<<""<<dEdxCut<<G4endl;
  return dEdxCut ;
}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::ComputeDEDX(const G4MaterialCutsCouple* matCC,
                                 const G4ParticleDefinition* p,
                                       G4double kineticEnergy,
                                       G4double cutEnergy)
{
  G4int iTkin,iPlace;
  size_t jMat;
  G4double massRatio  = fMass/p->GetPDGMass();
  G4double scaledTkin = kineticEnergy*massRatio;
  G4double charge     = p->GetPDGCharge();
  G4double charge2    = charge*charge, dEdx;

  for( jMat = 0 ;jMat < fMaterialCutsCoupleVector.size() ; ++jMat )
  {
    if( matCC == fMaterialCutsCoupleVector[jMat] ) break;
  }
  if(jMat == fMaterialCutsCoupleVector.size() && jMat > 0) jMat--;

  fPAIdEdxTable = fPAIdEdxBank[jMat];
  fdEdxVector = fdEdxTable[jMat];
  for(iTkin = 0 ; iTkin < fTotBin ; iTkin++)
  {
    if(scaledTkin < fParticleEnergyVector->GetLowEdgeEnergy(iTkin)) break ;    
  }
  iPlace = iTkin - 1;
  if(iPlace < 0) iPlace = 0;
  dEdx = charge2*( (*fdEdxVector)(iPlace) - GetdEdxCut(iPlace,cutEnergy) ) ;  

  if( dEdx < 0.) dEdx = 0.;
  return dEdx;
}

/////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::CrossSection( const G4MaterialCutsCouple* matCC,
                                   const G4ParticleDefinition* p,
                                         G4double kineticEnergy,
                                         G4double cutEnergy,
                                         G4double maxEnergy  ) 
{
  G4int iTkin,iPlace;
  size_t jMat;
  G4double tmax = min(MaxSecondaryEnergy(p, kineticEnergy), maxEnergy);
  G4double massRatio  = fMass/p->GetPDGMass();
  G4double scaledTkin = kineticEnergy*massRatio;
  G4double charge     = p->GetPDGCharge();
  G4double charge2    = charge*charge, cross, cross1, cross2;

  for( jMat = 0 ;jMat < fMaterialCutsCoupleVector.size() ; ++jMat )
  {
    if( matCC == fMaterialCutsCoupleVector[jMat] ) break;
  }
  if(jMat == fMaterialCutsCoupleVector.size() && jMat > 0) jMat--;

  fPAItransferTable = fPAIxscBank[jMat];

  for(iTkin = 0 ; iTkin < fTotBin ; iTkin++)
  {
    if(scaledTkin < fParticleEnergyVector->GetLowEdgeEnergy(iTkin)) break ;    
  }
  iPlace = iTkin - 1;
  if(iPlace < 0) iPlace = 0;

  // G4cout<<"iPlace = "<<iPlace<<"; tmax = "
  // <<tmax<<"; cutEnergy = "<<cutEnergy<<G4endl;  
  cross1 = GetdNdxCut(iPlace,tmax) ;
  // G4cout<<"cross1 = "<<cross1<<G4endl;  
  cross2 = GetdNdxCut(iPlace,cutEnergy) ;  
  // G4cout<<"cross2 = "<<cross2<<G4endl;  
  cross  = (cross2-cross1)*charge2;
  // G4cout<<"cross = "<<cross<<G4endl;  
  if( cross < DBL_MIN) cross = DBL_MIN;
  //  if( cross2 < DBL_MIN) cross2 = DBL_MIN;

  // return cross2;
  return cross;
}

///////////////////////////////////////////////////////////////////////////
//
// It is analog of PostStepDoIt in terms of secondary electron.
//

std::vector<G4DynamicParticle*>* 
G4PAIModel::SampleSecondaries( const G4MaterialCutsCouple* matCC,
                               const G4DynamicParticle* dp,
                                     G4double tmin,
                                     G4double maxEnergy)
{
  size_t jMat;
  for( jMat = 0 ;jMat < fMaterialCutsCoupleVector.size() ; ++jMat )
  {
    if( matCC == fMaterialCutsCoupleVector[jMat] ) break;
  }
  if(jMat == fMaterialCutsCoupleVector.size() && jMat > 0) jMat--;

  fPAItransferTable = fPAIxscBank[jMat];
  fdNdxCutVector    = fdNdxCutTable[jMat];

  G4double tmax = min(MaxSecondaryKinEnergy(dp), maxEnergy);
  if( tmin >= tmax && fVerbose)
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

  if( deltaTkin <= 0. && fVerbose ) 
  {
    G4cout<<"Tkin of secondary e- <= 0."<<G4endl;
    G4cout<<"G4PAIModel::SampleSecondary::deltaTkin = "<<deltaTkin<<G4endl;
    G4cout<<"G4PAIModel::SampleSecondary::deltaTkin = "<<deltaTkin<<G4endl;
    // deltaTkin = 10*eV;
    G4cout<<"Set G4PAIModel::SampleSecondary::deltaTkin = "<<deltaTkin<<G4endl;
  }
  if( deltaTkin <= 0.) return 0;

  if(deltaTkin > kineticEnergy && 
     particleMass != electron_mass_c2) deltaTkin = kineticEnergy;
  if (deltaTkin > 0.5*kineticEnergy && 
     (dp->GetDefinition()->GetParticleName() == "e-" || 
     dp->GetDefinition()->GetParticleName() == "e+")    ) deltaTkin = 0.5*kineticEnergy;

  G4double deltaTotalMomentum = sqrt(deltaTkin*(deltaTkin + 2. * electron_mass_c2 ));
  G4double totalMomentum      = sqrt(pSquare);
  G4double costheta           = deltaTkin*(totalEnergy + electron_mass_c2)
                                /(deltaTotalMomentum * totalMomentum);

  if( costheta >= 0.99999 ) costheta = 0.99999;
  G4double sintheta, sin2 = 1. - costheta*costheta;

  if( sin2 <= 0.) sintheta = 0.;
  else            sintheta = sqrt(sin2);

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

  std::vector<G4DynamicParticle*>* vdp = new std::vector<G4DynamicParticle*>;
  vdp->push_back(deltaRay);
  return vdp;
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

  for(iTkin=0;iTkin<fTotBin;iTkin++)
  {
    if(scaledTkin < fParticleEnergyVector->GetLowEdgeEnergy(iTkin))  break ;
  }
  iPlace = iTkin - 1 ;
  // G4cout<<"from search, iPlace = "<<iPlace<<G4endl ;
  if(iPlace < 0) iPlace = 0;
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

      for( iTransfer = 0;
 iTransfer < G4int((*fPAItransferTable)(iPlace)->GetVectorLength()); iTransfer++ )
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
  G4double x1, x2, y1, y2, energyTransfer ;

  if(iTransfer == 0)
  {
    energyTransfer = (*fPAItransferTable)(iPlace)->GetLowEdgeEnergy(iTransfer) ;
  }  
  else
  {
    if ( iTransfer >= G4int((*fPAItransferTable)(iPlace)->GetVectorLength()) )
    {
      iTransfer = (*fPAItransferTable)(iPlace)->GetVectorLength() - 1 ;
    }
    y1 = (*(*fPAItransferTable)(iPlace))(iTransfer-1) ;
    y2 = (*(*fPAItransferTable)(iPlace))(iTransfer) ;

    x1 = (*fPAItransferTable)(iPlace)->GetLowEdgeEnergy(iTransfer-1) ;
    x2 = (*fPAItransferTable)(iPlace)->GetLowEdgeEnergy(iTransfer) ;

    if ( x1 == x2 )    energyTransfer = x2 ;
    else
    {
      if ( y1 == y2  ) energyTransfer = x1 + (x2 - x1)*G4UniformRand() ;
      else
      {
        energyTransfer = x1 + (position - y1)*(x2 - x1)/(y2 - y1) ;
      }
    }
  }
  return energyTransfer ;
}

///////////////////////////////////////////////////////////////////////

G4double G4PAIModel::SampleFluctuations( const G4Material* material,
                                         const G4DynamicParticle* aParticle,
				               G4double&,
					       G4double& step,
                                               G4double&)
{
  size_t jMat;
  for( jMat = 0 ;jMat < fMaterialCutsCoupleVector.size() ; ++jMat )
  {
    if( material == fMaterialCutsCoupleVector[jMat]->GetMaterial() ) break;
  }
  if(jMat == fMaterialCutsCoupleVector.size() && jMat > 0) jMat--;

  fPAItransferTable = fPAIxscBank[jMat];
  fdNdxCutVector   = fdNdxCutTable[jMat];

  G4int iTkin, iTransfer, iPlace  ;
  G4long numOfCollisions=0;

  // G4cout<<"G4PAIModel::SampleFluctuations"<<G4endl ;
  // G4cout<<"in:  "<<fMaterialCutsCoupleVector[jMat]->GetMaterial()->GetName()<<G4endl ;

  G4double loss = 0.0, charge2 ;
  G4double stepSum = 0., stepDelta, lambda, omega; 
  G4double position, E1, E2, W1, W2, W, dNdxCut1, dNdxCut2, meanNumber;
  G4bool numb = true;
  G4double Tkin       = aParticle->GetKineticEnergy() ;
  G4double MassRatio  = fMass/aParticle->GetDefinition()->GetPDGMass() ;
  G4double charge     = aParticle->GetDefinition()->GetPDGCharge() ;
  charge2             = charge*charge ;
  G4double TkinScaled = Tkin*MassRatio ;

  for(iTkin=0;iTkin<fTotBin;iTkin++)
  {
    if(TkinScaled < fParticleEnergyVector->GetLowEdgeEnergy(iTkin))   break ;
  }
  iPlace = iTkin - 1 ; 
  if(iPlace < 0) iPlace = 0;
  //  G4cout<<"from search, iPlace = "<<iPlace<<G4endl ;
  dNdxCut1 = (*fdNdxCutVector)(iPlace) ;  
  //  G4cout<<"dNdxCut1 = "<<dNdxCut1<<G4endl ;


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
    //     G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl ;

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
      // G4cout<<"G4PAIModel::SampleFluctuations, omega = "<<omega/keV<<" keV; "<<"\t";
      loss += omega;
      numOfCollisions-- ;
    }
  }
  else
  {
    dNdxCut2 = (*fdNdxCutVector)(iPlace+1) ; 
    //  G4cout<<"dNdxCut2 = "<<dNdxCut2<<G4endl ;
 
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

      //  G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl ;

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
        // G4cout<<omega/keV<<"\t";
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

      // G4cout<<"(*(*fPAItransferTable)(iPlace))(0) = "<<
      //   (*(*fPAItransferTable)(iPlace))(0)<<G4endl ;
      // G4cout<<"(*(*fPAItransferTable)(iPlace+1))(0) = "<<
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

      //  G4cout<<"numOfCollisions = "<<numOfCollisions<<endl ;

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
  // G4cout<<"PAIModel AlongStepLoss = "<<loss/keV<<" keV, on step = "<<step/mm<<" mm"<<G4endl ; 
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


//
//
/////////////////////////////////////////////////






