//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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

const G4double G4PAIModel::fLowestKineticEnergy  = 10.0*keV ; // 100.0*MeV  ;
const G4double G4PAIModel::fHighestKineticEnergy = 100.*TeV ; // 1000.0*MeV ;
G4int G4PAIModel::fTotBin                        = 200  ;  // 50
G4int G4PAIModel::fMeanNumber                    = 20  ;  // 50 for dispersion

      // create physics vector and fill it

G4PhysicsLogVector*
G4PAIModel::fProtonEnergyVector = new G4PhysicsLogVector(fLowestKineticEnergy,
							   fHighestKineticEnergy,
							   fTotBin);

////////////////////////////////////////////////////////////////////////

G4PAIModel::G4PAIModel(const G4ParticleDefinition* p, const G4String& nam)
  : G4VEmModel(nam),G4VEmFluctuationModel(nam),
  fParticle(0),
  fHighKinEnergy(100.*TeV),
  fLowKinEnergy(2.0*MeV),
  fTwoln10(2.0*log(10.0)),
  fBg2lim(0.0169),
  fTaulim(8.4146e-3)
{
  if(p) SetParticle(p);
}

////////////////////////////////////////////////////////////////////////////

G4PAIModel::~G4PAIModel()
{
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

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::HighEnergyLimit(const G4ParticleDefinition* p)
{
  if(!fParticle) SetParticle(p);
  return fHighKinEnergy;
}

///////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::LowEnergyLimit( const G4ParticleDefinition* p )
{
  if(!fParticle) SetParticle(p);
  return fLowKinEnergy;
}

////////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::MinEnergyCut( const G4ParticleDefinition*,
                                   const G4MaterialCutsCouple* couple )
{
  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

////////////////////////////////////////////////////////////////////////////

G4bool G4PAIModel::IsInCharge( const G4ParticleDefinition* p )
{
  if(!fParticle) SetParticle(p);
  return (p->GetPDGCharge() != 0.0 );
}

////////////////////////////////////////////////////////////////////////////

void G4PAIModel::Initialise(const G4ParticleDefinition* p,
                                   const G4DataVector&)
{
  if(!fParticle) SetParticle(p);

  const G4ProductionCutsTable* theCoupleTable =
        G4ProductionCutsTable::GetProductionCutsTable();

  for(size_t iReg = 0; iReg < fPAIRegionVector.size();++iReg) // region loop
  {
    G4Region* curReg = fPAIRegionVector[iReg];

    // (*fPAIRegionVector[iRegion])

    std::vector<G4Material*>::const_iterator matIter = curReg->GetMaterialIterator();
    size_t jMat, numOfMat = curReg->GetNumberOfMaterials();

    //  for(size_t jMat = 0; jMat < curReg->GetNumberOfMaterials();++jMat){}

    for(jMat = 0 ; jMat < numOfMat; ++jMat) // region material loop
    {
      G4MaterialCutsCouple* matCouple = theCoupleTable->
      GetMaterialCutsCouple( *matIter, curReg->GetProductionCuts() );
      fMaterialCutsCoupleVector.push_back(matCouple);

      static const G4MaterialTable* theMaterialTable = 
                     G4Material::GetMaterialTable();
      G4int numberOfMat = G4Material::GetNumberOfMaterials();
      size_t iMatGlob;
      for(iMatGlob = 0 ; iMatGlob < numberOfMat ; iMatGlob++ )
      {
        if( *matIter == (*theMaterialTable)[iMatGlob]) break ;    
      }
      fMatIndex = iMatGlob;

      ComputeSandiaPhotoAbsCof();
      BuildPAIonisationTable();
      fPAIxscBank.push_back(fPAItransferBank);
      fPAIdEdxBank.push_back(fPAIdEdxTable);
      fdEdxTable.push_back(fdEdxVector);

      BuildLambdaVector(matCouple);
      fdNdxCutTable.push_back(fdNdxCutVector);
      fLambdaTable.push_back(fLambdaVector);

      // delete member used ???
      if(fdEdxVector)      delete fdEdxVector ;
      if ( fLambdaVector)  delete fLambdaVector;
      if ( fdNdxCutVector) delete fdNdxCutVector;
      if( fPAItransferBank )    
      {
        fPAItransferBank->clearAndDestroy();
        delete fPAItransferBank ;
      }
      for(G4int i=0;i<fSandiaIntervalNumber;i++)
      {
        delete[] fSandiaPhotoAbsCof[i];
      }
      delete[] fSandiaPhotoAbsCof;

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

  for(i=0;i<numberOfElements;i++)  thisMaterialZ[i] = 
    (G4int)(*theMaterialTable)[fMatIndex]->GetElement(i)->GetZ() ;
   
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
  delete[] thisMaterialZ ;
}

////////////////////////////////////////////////////////////////////////////
//
// Build tables for the ionization energy loss
//  the tables are built for MATERIALS
//                           *********

void
G4PAIModel::BuildPAIonisationTable()
{
  G4double LowEdgeEnergy , ionloss ;
  G4double massRatio, tau, Tmax, Tmin, Tkin, deltaLow, gamma, bg2 ;
   
  if( fPAItransferBank )    
  {
     fPAItransferBank->clearAndDestroy() ;
     delete fPAItransferBank ;
  }
  fPAItransferBank = new G4PhysicsTable(fTotBin);
        
  if(fdEdxVector) delete fdEdxVector ;
  fdEdxVector = new G4PhysicsLogVector( fLowestKineticEnergy, 
					 fHighestKineticEnergy,
					 fTotBin               ) ;   
  Tmin     = fSandiaPhotoAbsCof[0][0] ;      // low energy Sandia interval
  deltaLow = 0.5*eV ;         

  for (G4int i = 0 ; i < fTotBin ; i++)  //The loop for the kinetic energy
  {
    LowEdgeEnergy = fProtonEnergyVector->GetLowEdgeEnergy(i) ;
    tau = LowEdgeEnergy/proton_mass_c2 ;
    //    if(tau < 0.01)  tau = 0.01 ;
    gamma = tau +1. ;
    // G4cout<<"gamma = "<<gamma<<endl ;
    bg2 = tau*(tau + 2. ) ;
    massRatio = electron_mass_c2/proton_mass_c2 ;
    Tmax = 2.*electron_mass_c2*bg2/(1.+2.*gamma*massRatio+massRatio*massRatio) ;
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

    fPAItransferBank->insertAt(i,transferVector) ;
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

  size_t numOfCouples = theCoupleTable->GetTableSize(), jMatCC;
  for (jMatCC = 0 ; jMatCC < numOfCouples ; jMatCC++ )
  {
    if( matCutsCouple == theCoupleTable->GetMaterialCutsCouple(jMatCC) ) break;
  } 


  const std::vector<G4double>*  deltaCutInKineticEnergy = theCoupleTable->
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

  G4cout<<"PAIModel DeltaCutInKineticEnergyNow = "
        <<deltaCutInKineticEnergyNow/keV<<" keV"<<G4endl;

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
  // G4cout<<"size = "<<G4int((*fPAItransferBank)(iPlace)->GetVectorLength())
  //           <<G4endl;  
  for( iTransfer = 0 ; 
       iTransfer < G4int((*fPAItransferBank)(iPlace)->GetVectorLength()) ; 
       iTransfer++)
  {
    if(transferCut <= (*fPAItransferBank)(iPlace)->GetLowEdgeEnergy(iTransfer))
    {
      break ;
    }
  }  
  if ( iTransfer >= G4int((*fPAItransferBank)(iPlace)->GetVectorLength()) )
  {
      iTransfer = (*fPAItransferBank)(iPlace)->GetVectorLength() - 1 ;
  }
  y1 = (*(*fPAItransferBank)(iPlace))(iTransfer-1) ;
  y2 = (*(*fPAItransferBank)(iPlace))(iTransfer) ;
  // G4cout<<"y1 = "<<y1<<"; "<<"y2 = "<<y2<<G4endl;
  x1 = (*fPAItransferBank)(iPlace)->GetLowEdgeEnergy(iTransfer-1) ;
  x2 = (*fPAItransferBank)(iPlace)->GetLowEdgeEnergy(iTransfer) ;
  // G4cout<<"x1 = "<<x1<<"; "<<"x2 = "<<x2<<G4endl;

  if ( y1 == y2 )    dNdxCut = y2 ;
  else
  {
    //  if ( x1 == x2  ) dNdxCut = y1 + (y2 - y1)*G4UniformRand() ;
    if ( abs(x1-x2) <= eV  ) dNdxCut = y1 + (y2 - y1)*G4UniformRand() ;
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
    if ( abs(x1-x2) <= eV  ) dEdxCut = y1 + (y2 - y1)*G4UniformRand() ;
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
  G4double scaledTkin = kineticEnergy*p->GetPDGMass()/proton_mass_c2;
  G4double charge     = p->GetPDGCharge();
  G4double charge2    = charge*charge, dEdx;

  for( jMat = 0 ;jMat < fMaterialCutsCoupleVector.size() ; ++jMat )
  {
    if( matCC == fMaterialCutsCoupleVector[jMat] ) break;
  }
  fPAIdEdxTable = fPAIdEdxBank[jMat];
  fdEdxVector = fdEdxTable[jMat];
  for(iTkin = 0 ; iTkin < fTotBin ; iTkin++)
  {
    if(scaledTkin < fProtonEnergyVector->GetLowEdgeEnergy(iTkin)) break ;    
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
  G4double tmax = std::min(MaxSecondaryEnergy(p, kineticEnergy), maxEnergy);
  G4double scaledTkin = kineticEnergy*p->GetPDGMass()/proton_mass_c2;
  G4double charge     = p->GetPDGCharge();
  G4double charge2    = charge*charge, cross, cross1, cross2;

  for( jMat = 0 ;jMat < fMaterialCutsCoupleVector.size() ; ++jMat )
  {
    if( matCC == fMaterialCutsCoupleVector[jMat] ) break;
  }
  fPAItransferBank = fPAIxscBank[jMat];

  for(iTkin = 0 ; iTkin < fTotBin ; iTkin++)
  {
    if(scaledTkin < fProtonEnergyVector->GetLowEdgeEnergy(iTkin)) break ;    
  }
  iPlace = iTkin - 1;
  if(iPlace < 0) iPlace = 0;
  cross1 = GetdNdxCut(iPlace,tmax) ;  
  cross2 = GetdNdxCut(iPlace,cutEnergy) ;  
  cross  = (cross2-cross1)*charge2;
  if( cross < 0.) cross = 0.;
  return cross;
}

///////////////////////////////////////////////////////////////////////////
//
// It is analog of PostStepDoIt in terms of secondary electron.
//

G4DynamicParticle* 
G4PAIModel::SampleSecondary( const G4MaterialCutsCouple* matCC,
                             const G4DynamicParticle* dp,
                                   G4double tmin,
                                   G4double maxEnergy)
{
  size_t jMat;
  for( jMat = 0 ;jMat < fMaterialCutsCoupleVector.size() ; ++jMat )
  {
    if( matCC == fMaterialCutsCoupleVector[jMat] ) break;
  }
  fPAItransferBank = fPAIxscBank[jMat];
  fdNdxCutVector     = fdNdxCutTable[jMat];

  G4double tmax = std::min(MaxSecondaryEnergy(dp), maxEnergy);
  if( tmin >= tmax ) return 0;
  G4ThreeVector momentum = dp->GetMomentumDirection();
  G4double particleMass  = dp->GetMass();
  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double scaledTkin    = kineticEnergy*particleMass/proton_mass_c2;
  G4double totalEnergy   = kineticEnergy + particleMass;
  G4double pSquare       = kineticEnergy*(totalEnergy+particleMass);
 
  G4double deltaTkin     = GetPostStepTransfer(scaledTkin);
  if( deltaTkin <= 0. ) return 0;
  G4double deltaTotalMomentum = sqrt(deltaTkin*(deltaTkin + 2. * electron_mass_c2 ));
  G4double totalMomentum      = sqrt(pSquare);
  G4double costheta           = deltaTkin*(totalEnergy + electron_mass_c2)
                                /(deltaTotalMomentum * totalMomentum);
  if (costheta < 0.) costheta = 0.;
  if (costheta > +1.) costheta = +1.;

  //  direction of the delta electron
  
  G4double phi = twopi*G4UniformRand(); 
  G4double sintheta = sqrt((1.+costheta)*(1.-costheta));
  G4double dirx = sintheta*cos(phi), diry = sintheta*sin(phi), dirz = costheta;

  G4ThreeVector deltaDirection(dirx,diry,dirz);
  deltaDirection.rotateUz(momentum);

  // create G4DynamicParticle object for delta ray
 
  G4DynamicParticle* deltaRay = new G4DynamicParticle;
  deltaRay->SetDefinition(G4Electron::Electron());
  deltaRay->SetKineticEnergy( deltaTkin );
  deltaRay->SetMomentumDirection(deltaDirection); 

  return deltaRay;
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
    if(scaledTkin < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))  break ;
  }
  iPlace = iTkin - 1 ;
  dNdxCut1 = (*fdNdxCutVector)(iPlace) ;  

  //  G4cout<<"iPlace = "<<iPlace<<endl ;

  if(iTkin == fTotBin) // Fermi plato, try from left
  {
      position = dNdxCut1*G4UniformRand() ;

      for( iTransfer = 0;
 iTransfer < G4int((*fPAItransferBank)(iPlace)->GetVectorLength()); iTransfer++ )
      {
        if(position >= (*(*fPAItransferBank)(iPlace))(iTransfer)) break ;
      }
      transfer = GetEnergyTransfer(iPlace,position,iTransfer);
  }
  else
  {
    dNdxCut2 = (*fdNdxCutVector)(iPlace+1) ;  
    if(iTkin == 0) // Tkin is too small, trying from right only
    {
      position = dNdxCut2*G4UniformRand() ;

      for( iTransfer = 0;
  iTransfer < G4int((*fPAItransferBank)(iPlace+1)->GetVectorLength()); iTransfer++ )
      {
        if(position >= (*(*fPAItransferBank)(iPlace+1))(iTransfer)) break ;
      }
      transfer = GetEnergyTransfer(iPlace+1,position,iTransfer);
    } 
    else // general case: Tkin between two vectors of the material
    {
      E1 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1) ; 
      E2 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin)     ;
      W  = 1.0/(E2 - E1) ;
      W1 = (E2 - scaledTkin)*W ;
      W2 = (scaledTkin - E1)*W ;

      position = ( dNdxCut1*W1 + dNdxCut2*W2 )*G4UniformRand() ;

        // G4cout<<position<<"\t" ;

      for( iTransfer = 0;
 iTransfer < G4int((*fPAItransferBank)(iPlace)->GetVectorLength()); iTransfer++ )
      {
          if( position >=
          ( (*(*fPAItransferBank)(iPlace))(iTransfer)*W1 + 
            (*(*fPAItransferBank)(iPlace+1))(iTransfer)*W2) ) break ;
      }
      transfer = GetEnergyTransfer(iPlace,position,iTransfer);
    }
  } 
  // G4cout<<"PAImodel PostStepTransfer = "<<transfer/keV<<" keV"<<endl ; 
  if(transfer < 0.0 ) transfer = 0.0 ;
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
    energyTransfer = (*fPAItransferBank)(iPlace)->GetLowEdgeEnergy(iTransfer) ;
  }  
  else
  {
    if ( iTransfer >= G4int((*fPAItransferBank)(iPlace)->GetVectorLength()) )
    {
      iTransfer = (*fPAItransferBank)(iPlace)->GetVectorLength() - 1 ;
    }
    y1 = (*(*fPAItransferBank)(iPlace))(iTransfer-1) ;
    y2 = (*(*fPAItransferBank)(iPlace))(iTransfer) ;

    x1 = (*fPAItransferBank)(iPlace)->GetLowEdgeEnergy(iTransfer-1) ;
    x2 = (*fPAItransferBank)(iPlace)->GetLowEdgeEnergy(iTransfer) ;

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


////////////////////////////////////////////////////////////////////////////

std::vector<G4DynamicParticle*>* 
G4PAIModel::SampleSecondaries( const G4MaterialCutsCouple* couple,
                               const G4DynamicParticle* dp,
                                     G4double tmin,
                                     G4double maxEnergy)
{
  std::vector<G4DynamicParticle*>* vdp = new std::vector<G4DynamicParticle*>;
  G4DynamicParticle* delta             = SampleSecondary(couple, dp, tmin, maxEnergy);
  vdp->push_back(delta);
  return vdp;
}

///////////////////////////////////////////////////////////////////////

G4double G4PAIModel::SampleFluctuations( const G4Material* material,
                                         const G4DynamicParticle* aParticle,
				               G4double& tmax,
					       G4double& step,
                                               G4double& meanLoss     )
{
  size_t jMat;
  for( jMat = 0 ;jMat < fMaterialCutsCoupleVector.size() ; ++jMat )
  {
    if( material == fMaterialCutsCoupleVector[jMat]->GetMaterial() ) break;
  }
  fPAItransferBank = fPAIxscBank[jMat];
  fdNdxCutVector   = fdNdxCutTable[jMat];

  G4int iTkin, iTransfer, iPlace  ;
  G4long numOfCollisions;

  // G4cout<<"G4PAIModel::SampleFluctuations"<<G4endl ;

  G4double loss = 0.0, charge2 ;
 
  G4double position, E1, E2, W1, W2, W, dNdxCut1, dNdxCut2, meanNumber;

  G4double Tkin       = aParticle->GetKineticEnergy() ;
  G4double MassRatio  = proton_mass_c2/aParticle->GetDefinition()->GetPDGMass() ;
  G4double charge     = aParticle->GetDefinition()->GetPDGCharge() ;
  charge2             = charge*charge ;
  G4double TkinScaled = Tkin*MassRatio ;

  for(iTkin=0;iTkin<fTotBin;iTkin++)
  {
    if(TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))   break ;
  }
  iPlace = iTkin - 1 ; 
  dNdxCut1 = (*fdNdxCutVector)(iPlace) ;  

  //  G4cout<<"iPlace = "<<iPlace<<endl ;

  if(iTkin == fTotBin) // Fermi plato, try from left
  {
    meanNumber =((*(*fPAItransferBank)(iPlace))(0)-dNdxCut1)*step*charge2;
    if(meanNumber < 0.) meanNumber = 0. ;
    numOfCollisions = RandPoisson::shoot(meanNumber) ;
    
    //     G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl ;

    while(numOfCollisions)
    {
      position = dNdxCut1+
                 ((*(*fPAItransferBank)(iPlace))(0)-dNdxCut1)*G4UniformRand() ;

      for( iTransfer = 0;
   iTransfer < G4int((*fPAItransferBank)(iPlace)->GetVectorLength()); iTransfer++ )
      {
        if(position >= (*(*fPAItransferBank)(iPlace))(iTransfer)) break ;
      }
      loss += GetEnergyTransfer(iPlace,position,iTransfer);
      numOfCollisions-- ;
    }
  }
  else
  {
    dNdxCut2 = (*fdNdxCutVector)(iPlace+1) ; 
 
    if(iTkin == 0) // Tkin is too small, trying from right only
    {
      meanNumber =((*(*fPAItransferBank)(iPlace+1))(0)-dNdxCut2)*step*charge2;
      if( meanNumber < 0. ) meanNumber = 0. ;
      numOfCollisions = RandPoisson::shoot(meanNumber) ;

      //  G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl ;

      while(numOfCollisions)
      {
        position = dNdxCut2+
                   ((*(*fPAItransferBank)(iPlace+1))(0)-dNdxCut2)*G4UniformRand();
   
        for( iTransfer = 0;
   iTransfer < G4int((*fPAItransferBank)(iPlace+1)->GetVectorLength()); iTransfer++ )
        {
          if(position >= (*(*fPAItransferBank)(iPlace+1))(iTransfer)) break ;
        }
        loss += GetEnergyTransfer(iPlace+1,position,iTransfer);
        numOfCollisions-- ;
      }
    } 
    else // general case: Tkin between two vectors of the material
    {
      E1 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1) ; 
      E2 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin)     ;
       W = 1.0/(E2 - E1) ;
      W1 = (E2 - TkinScaled)*W ;
      W2 = (TkinScaled - E1)*W ;

      // G4cout<<"(*(*fPAItransferBank)(iPlace))(0) = "<<
      //   (*(*fPAItransferBank)(iPlace))(0)<<G4endl ;
      // G4cout<<"(*(*fPAItransferBank)(iPlace+1))(0) = "<<
      //     (*(*fPAItransferBank)(iPlace+1))(0)<<G4endl ;

      meanNumber=( ((*(*fPAItransferBank)(iPlace))(0)-dNdxCut1)*W1 + 
		   ((*(*fPAItransferBank)(iPlace+1))(0)-dNdxCut2)*W2 )*step*charge2;
      if(meanNumber<0.0) meanNumber = 0.0;
      numOfCollisions = RandPoisson::shoot(meanNumber) ;

      //  G4cout<<"numOfCollisions = "<<numOfCollisions<<endl ;

      while(numOfCollisions)
      {
        position =( (dNdxCut1+
                  ((*(*fPAItransferBank)(iPlace  ))(0)-dNdxCut1))*W1 + 
                    (dNdxCut2+
                  ((*(*fPAItransferBank)(iPlace+1))(0)-dNdxCut2))*W2 )*G4UniformRand();

        // G4cout<<position<<"\t" ;

        for( iTransfer = 0;
    iTransfer < G4int((*fPAItransferBank)(iPlace)->GetVectorLength()); iTransfer++ )
        {
          if( position >=
          ( (*(*fPAItransferBank)(iPlace))(iTransfer)*W1 + 
            (*(*fPAItransferBank)(iPlace+1))(iTransfer)*W2) )
          {
	      break ;
	  }
        }
	// loss += (*fPAItransferBank)(iPlace)->GetLowEdgeEnergy(iTransfer) ; 
        loss += GetEnergyTransfer(iPlace,position,iTransfer);
        numOfCollisions-- ;    
      }
    }
  } 
  // G4cout<<"PAIModel AlongStepLoss = "<<loss/keV<<" keV"<<endl ; 

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






