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

// ------------------------------------------------------------
// G4RDHadronIonisation
//
//
// Author: Maria Grazia Pia (MariaGrazia.Pia@ge.infn.it)
//
// 08 Sep 2008 - MGP - Created (initially based on G4hLowEnergyIonisation) 
//                     Added PIXE capabilities
//                     Partial clean-up of the implementation (more needed)
//                     Calculation of MicroscopicCrossSection delegated to specialised cla// Documentation available in:
// M.G. Pia et al., PIXE Simulation With Geant4,
// IEEE Trans. Nucl. Sci., vol. 56, no. 6, pp. 3614-3649, Dec. 2009.

//
// ------------------------------------------------------------
 
#include "G4hImpactIonisation.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Poisson.hh"
#include "G4UnitsTable.hh"
#include "G4EnergyLossTables.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4AtomicDeexcitation.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4PixeCrossSectionHandler.hh"
#include "G4IInterpolator.hh"
#include "G4LogLogInterpolator.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"               
#include "G4ProcessManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PhysicsLogVector.hh"       
#include "G4PhysicsLinearVector.hh"    

#include "G4VLowEnergyModel.hh"
#include "G4hNuclearStoppingModel.hh"  
#include "G4hBetheBlochModel.hh"       
#include "G4hParametrisedLossModel.hh" 
#include "G4QAOLowEnergyLoss.hh"       
#include "G4hIonEffChargeSquare.hh"    
#include "G4IonChuFluctuationModel.hh" 
#include "G4IonYangFluctuationModel.hh"

#include "G4MaterialCutsCouple.hh"
#include "G4Track.hh"
#include "G4Step.hh"

G4hImpactIonisation::G4hImpactIonisation(const G4String& processName)
  : G4hRDEnergyLoss(processName),
    betheBlochModel(0),
    protonModel(0),
    antiprotonModel(0),
    theIonEffChargeModel(0),
    theNuclearStoppingModel(0),
    theIonChuFluctuationModel(0),
    theIonYangFluctuationModel(0),
    protonTable("ICRU_R49p"),
    antiprotonTable("ICRU_R49p"),
    theNuclearTable("ICRU_R49"),
    nStopping(true),
    theBarkas(true),
    theMeanFreePathTable(0),
    paramStepLimit (0.005),
    pixeCrossSectionHandler(0)
{ 
  InitializeMe();
}



void G4hImpactIonisation::InitializeMe()
{
  LowestKineticEnergy  = 10.0*eV ;
  HighestKineticEnergy = 100.0*GeV ;
  MinKineticEnergy     = 10.0*eV ; 
  TotBin               = 360 ;
  protonLowEnergy      = 1.*keV ;
  protonHighEnergy     = 100.*MeV ;
  antiprotonLowEnergy  = 25.*keV ;
  antiprotonHighEnergy = 2.*MeV ;
  minGammaEnergy       = 100 * eV;
  minElectronEnergy    = 250.* eV;
  verboseLevel         = 0;
 
  // Min and max energy of incident particle for the calculation of shell cross sections
  // for PIXE generation
  eMinPixe = 1.* keV;
  eMaxPixe = 200. * MeV;
  
  G4String defaultPixeModel("ecpssr"); 
  modelK = defaultPixeModel;
  modelL = defaultPixeModel;
  modelM = defaultPixeModel;
}


G4hImpactIonisation::~G4hImpactIonisation()
{
  if (theMeanFreePathTable) 
    {
      theMeanFreePathTable->clearAndDestroy();
      delete theMeanFreePathTable;
    }

  if (betheBlochModel) delete betheBlochModel;
  if (protonModel) delete protonModel;
  if (antiprotonModel) delete antiprotonModel;
  if (theNuclearStoppingModel) delete theNuclearStoppingModel;
  if (theIonEffChargeModel) delete theIonEffChargeModel;
  if (theIonChuFluctuationModel) delete theIonChuFluctuationModel;
  if (theIonYangFluctuationModel) delete theIonYangFluctuationModel;

  delete pixeCrossSectionHandler;

  // ---- MGP ---- The following is to be checked
  //  if (shellVacancy) delete shellVacancy;

  cutForDelta.clear();
}

// --------------------------------------------------------------------
void G4hImpactIonisation::SetElectronicStoppingPowerModel(const G4ParticleDefinition* particle,
							  const G4String& dedxTable)
// This method defines the ionisation parametrisation method via its name 
{
  if (particle->GetPDGCharge() > 0 ) 
    {
      // Positive charge
      SetProtonElectronicStoppingPowerModel(dedxTable) ;
    } 
  else 
    {
      // Antiprotons
      SetAntiProtonElectronicStoppingPowerModel(dedxTable) ;
    }
}


// --------------------------------------------------------------------
void G4hImpactIonisation::InitializeParametrisation() 

{
  // Define models for parametrisation of electronic energy losses
  betheBlochModel = new G4hBetheBlochModel("Bethe-Bloch") ;
  protonModel = new G4hParametrisedLossModel(protonTable) ;
  protonHighEnergy = std::min(protonHighEnergy,protonModel->HighEnergyLimit(0, 0));
  antiprotonModel = new G4QAOLowEnergyLoss(antiprotonTable) ;
  theNuclearStoppingModel = new G4hNuclearStoppingModel(theNuclearTable) ;
  theIonEffChargeModel = new G4hIonEffChargeSquare("Ziegler1988") ;
  theIonChuFluctuationModel = new G4IonChuFluctuationModel("Chu") ;
  theIonYangFluctuationModel = new G4IonYangFluctuationModel("Yang") ;
}


// --------------------------------------------------------------------
void G4hImpactIonisation::BuildPhysicsTable(const G4ParticleDefinition& particleDef)

//  just call BuildLossTable+BuildLambdaTable
{

  // Verbose print-out
  if(verboseLevel > 0) 
    {
      G4cout << "G4hImpactIonisation::BuildPhysicsTable for "
	     << particleDef.GetParticleName()
	     << " mass(MeV)= " << particleDef.GetPDGMass()/MeV
	     << " charge= " << particleDef.GetPDGCharge()/eplus
	     << " type= " << particleDef.GetParticleType()
	     << G4endl;
      
      if(verboseLevel > 1) 
	{
	  G4ProcessVector* pv = particleDef.GetProcessManager()->GetProcessList();
	  
	  G4cout << " 0: " << (*pv)[0]->GetProcessName() << " " << (*pv)[0]
		 << " 1: " << (*pv)[1]->GetProcessName() << " " << (*pv)[1]
	    //        << " 2: " << (*pv)[2]->GetProcessName() << " " << (*pv)[2]
		 << G4endl;
	  G4cout << "ionModel= " << theIonEffChargeModel
		 << " MFPtable= " << theMeanFreePathTable
		 << " iniMass= " << initialMass
		 << G4endl;
	}
    }
  // End of verbose print-out

  if (particleDef.GetParticleType() == "nucleus" &&
      particleDef.GetParticleName() != "GenericIon" &&
      particleDef.GetParticleSubType() == "generic")
    {

      G4EnergyLossTables::Register(&particleDef,
				   theDEDXpTable,
				   theRangepTable,
				   theInverseRangepTable,
				   theLabTimepTable,
				   theProperTimepTable,
				   LowestKineticEnergy, HighestKineticEnergy,
				   proton_mass_c2/particleDef.GetPDGMass(),
				   TotBin);

      return;
    }

  if( !CutsWhereModified() && theLossTable) return;

  InitializeParametrisation() ;
  G4Proton* proton = G4Proton::Proton();
  G4AntiProton* antiproton = G4AntiProton::AntiProton();

  charge = particleDef.GetPDGCharge() / eplus;
  chargeSquare = charge*charge ;

  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  cutForDelta.clear();
  cutForGamma.clear();

  for (G4int j=0; j<numOfCouples; ++j) {

    // get material parameters needed for the energy loss calculation
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);
    const G4Material* material= couple->GetMaterial();

    // the cut cannot be below lowest limit
    G4double tCut = (*(theCoupleTable->GetEnergyCutsVector(1)))[j];
    if(tCut > HighestKineticEnergy) tCut = HighestKineticEnergy;

    G4double excEnergy = material->GetIonisation()->GetMeanExcitationEnergy();

    tCut = std::max(tCut,excEnergy);
    cutForDelta.push_back(tCut);

    // the cut cannot be below lowest limit
    tCut = (*(theCoupleTable->GetEnergyCutsVector(0)))[j];
    if(tCut > HighestKineticEnergy) tCut = HighestKineticEnergy;
    tCut = std::max(tCut,minGammaEnergy);
    cutForGamma.push_back(tCut);
  }

  if(verboseLevel > 0) {
    G4cout << "Cuts are defined " << G4endl;
  }

  if(0.0 < charge)
    {
      {
        BuildLossTable(*proton) ;

	//      The following vector has a fixed dimension (see src/G4hImpactLoss.cc for more details)        
	//      It happended in the past that caused memory corruption errors. The problem is still pending, even if temporary solved
	//        G4cout << "[NOTE]: __LINE__=" << __LINE__ << ", particleDef=" << particleDef.GetParticleName() << ", proton=" << proton << ", theLossTable=" << theLossTable << ", CounterOfpProcess=" << CounterOfpProcess << G4endl;
        
        RecorderOfpProcess[CounterOfpProcess] = theLossTable ;
        CounterOfpProcess++;
      }
    } else {
    {
      BuildLossTable(*antiproton) ;
        
      //      The following vector has a fixed dimension (see src/G4hImpactLoss.cc for more details)        
      //      It happended in the past that caused memory corruption errors. The problem is still pending, even if temporary solved
      //        G4cout << "[NOTE]: __LINE__=" << __LINE__ << ", particleDef=" << particleDef.GetParticleName() << ", antiproton=" << antiproton << ", theLossTable=" << theLossTable << ", CounterOfpbarProcess=" << CounterOfpbarProcess << G4endl;
        
      RecorderOfpbarProcess[CounterOfpbarProcess] = theLossTable ;
      CounterOfpbarProcess++;
    }
  }

  if(verboseLevel > 0) {
    G4cout << "G4hImpactIonisation::BuildPhysicsTable: "
           << "Loss table is built "
      //	   << theLossTable
	   << G4endl;
  }

  BuildLambdaTable(particleDef) ;
  //  BuildDataForFluorescence(particleDef);

  if(verboseLevel > 1) {
    G4cout << (*theMeanFreePathTable) << G4endl;
  }

  if(verboseLevel > 0) {
    G4cout << "G4hImpactIonisation::BuildPhysicsTable: "
           << "DEDX table will be built "
      //	   << theDEDXpTable << " " << theDEDXpbarTable
      //	   << " " << theRangepTable << " " << theRangepbarTable
	   << G4endl;
  }

  BuildDEDXTable(particleDef) ;

  if(verboseLevel > 1) {
    G4cout << (*theDEDXpTable) << G4endl;
  }

  if((&particleDef == proton) ||  (&particleDef == antiproton)) PrintInfoDefinition() ;

  if(verboseLevel > 0) {
    G4cout << "G4hImpactIonisation::BuildPhysicsTable: end for "
           << particleDef.GetParticleName() << G4endl;
  }
}





// --------------------------------------------------------------------
void G4hImpactIonisation::BuildLossTable(const G4ParticleDefinition& particleDef)
{
  // Initialisation
  G4double lowEdgeEnergy , ionloss, ionlossBB, paramB ;
  //G4double lowEnergy, highEnergy;
  G4double highEnergy;
  G4Proton* proton = G4Proton::Proton();

  if (particleDef == *proton) 
    {
      //lowEnergy = protonLowEnergy ;
      highEnergy = protonHighEnergy ;
      charge = 1. ;
    } 
  else 
    {
      //lowEnergy = antiprotonLowEnergy ;
      highEnergy = antiprotonHighEnergy ;
      charge = -1. ;
    }
  chargeSquare = 1. ;

  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  if ( theLossTable) 
    {
      theLossTable->clearAndDestroy();
      delete theLossTable;
    }

  theLossTable = new G4PhysicsTable(numOfCouples);

  //  loop for materials
  for (G4int j=0; j<numOfCouples; ++j) {

    // create physics vector and fill it
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(LowestKineticEnergy,
                                                         HighestKineticEnergy,
                                                         TotBin);

    // get material parameters needed for the energy loss calculation
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);
    const G4Material* material= couple->GetMaterial();

    if ( charge > 0.0 ) {
      ionloss = ProtonParametrisedDEDX(couple,highEnergy) ;
    } else {
      ionloss = AntiProtonParametrisedDEDX(couple,highEnergy) ;
    }

    ionlossBB = betheBlochModel->TheValue(&particleDef,material,highEnergy) ;
    ionlossBB -= DeltaRaysEnergy(couple,highEnergy,proton_mass_c2) ;


    paramB =  ionloss/ionlossBB - 1.0 ;

    // now comes the loop for the kinetic energy values
    for (G4int i = 0 ; i < TotBin ; i++) {
      lowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;

      // low energy part for this material, parametrised energy loss formulae
      if ( lowEdgeEnergy < highEnergy ) {

        if ( charge > 0.0 ) {
          ionloss = ProtonParametrisedDEDX(couple,lowEdgeEnergy) ;
	} else {
          ionloss = AntiProtonParametrisedDEDX(couple,lowEdgeEnergy) ;
	}

      } else {

        // high energy part for this material, Bethe-Bloch formula
        ionloss = betheBlochModel->TheValue(proton,material,
					    lowEdgeEnergy) ;

        ionloss -= DeltaRaysEnergy(couple,lowEdgeEnergy,proton_mass_c2) ;

	ionloss *= (1.0 + paramB*highEnergy/lowEdgeEnergy) ;
      }

      // now put the loss into the vector
      if(verboseLevel > 1) {
        G4cout << "E(MeV)= " << lowEdgeEnergy/MeV
               << "  dE/dx(MeV/mm)= " << ionloss*mm/MeV
               << " in " << material->GetName() << G4endl;
      }
      aVector->PutValue(i,ionloss) ;
    }
    // Insert vector for this material into the table
    theLossTable->insert(aVector) ;
  }
}



// --------------------------------------------------------------------
void G4hImpactIonisation::BuildLambdaTable(const G4ParticleDefinition& particleDef)

{
  // Build mean free path tables for the delta ray production process
  //     tables are built for MATERIALS

  if(verboseLevel > 1) {
    G4cout << "G4hImpactIonisation::BuildLambdaTable for "
           << particleDef.GetParticleName() << " is started" << G4endl;
  }


  G4double lowEdgeEnergy, value;
  charge = particleDef.GetPDGCharge()/eplus ;
  chargeSquare = charge*charge ;
  initialMass = particleDef.GetPDGMass();

  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();


  if (theMeanFreePathTable) {
    theMeanFreePathTable->clearAndDestroy();
    delete theMeanFreePathTable;
  }

  theMeanFreePathTable = new G4PhysicsTable(numOfCouples);

  // loop for materials

  for (G4int j=0 ; j < numOfCouples; ++j) {

    //create physics vector then fill it ....
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(LowestKineticEnergy,
                                                         HighestKineticEnergy,
                                                         TotBin);

    // compute the (macroscopic) cross section first
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);
    const G4Material* material= couple->GetMaterial();

    const G4ElementVector* theElementVector =  material->GetElementVector() ;
    const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();
    const G4int numberOfElements = (G4int)material->GetNumberOfElements() ;

    // get the electron kinetic energy cut for the actual material,
    //  it will be used in ComputeMicroscopicCrossSection
    // ( it is the SAME for ALL the ELEMENTS in THIS MATERIAL )
    //   ------------------------------------------------------

    G4double deltaCut = cutForDelta[j];

    for ( G4int i = 0 ; i < TotBin ; i++ ) {
      lowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
      G4double sigma = 0.0 ;
      G4int Z;
      
      for (G4int iel=0; iel<numberOfElements; iel++ ) 
	{
	  Z = (G4int) (*theElementVector)[iel]->GetZ();
	  // ---- MGP --- Corrected duplicated cross section calculation here from
	  // G4hLowEnergyIonisation original
	  G4double microCross = MicroscopicCrossSection( particleDef,
							 lowEdgeEnergy,
							 Z,
							 deltaCut ) ;	
	  //totalCrossSectionMap [Z] = microCross;
	  sigma += theAtomicNumDensityVector[iel] * microCross ; 
	}

      // mean free path = 1./macroscopic cross section

      value = sigma<=0 ? DBL_MAX : 1./sigma ;

      aVector->PutValue(i, value) ;
    }

    theMeanFreePathTable->insert(aVector);
  }

}


// --------------------------------------------------------------------
G4double G4hImpactIonisation::MicroscopicCrossSection(const G4ParticleDefinition& particleDef,
						      G4double kineticEnergy,
						      G4double atomicNumber,
						      G4double deltaCutInEnergy) const
{
  //******************************************************************
    // cross section formula is OK for spin=0, 1/2, 1 only !
    // *****************************************************************

    // Calculates the microscopic cross section in GEANT4 internal units
    // Formula documented in Geant4 Phys. Ref. Manual
    // ( it is called for elements, AtomicNumber = z )

    G4double totalCrossSection = 0.;

  // Particle mass and energy
  G4double particleMass = initialMass;
  G4double energy = kineticEnergy + particleMass;

  // Some kinematics
  G4double gamma = energy / particleMass;
  G4double beta2 = 1. - 1. / (gamma * gamma);
  G4double var = electron_mass_c2 / particleMass;
  G4double tMax   = 2. * electron_mass_c2 * (gamma*gamma - 1.) / (1. + 2.* gamma*var + var*var);

  // Calculate the total cross section

  if ( tMax > deltaCutInEnergy ) 
    {
      var = deltaCutInEnergy / tMax;
      totalCrossSection = (1. - var * (1. - beta2 * std::log(var))) / deltaCutInEnergy ;
      
      G4double spin = particleDef.GetPDGSpin() ;
      
      // +term for spin=1/2 particle
      if (spin == 0.5) 
	{
	  totalCrossSection +=  0.5 * (tMax - deltaCutInEnergy) / (energy*energy);
	}
      // +term for spin=1 particle
      else if (spin > 0.9 )
	{
	  totalCrossSection += -std::log(var) / 
	    (3. * deltaCutInEnergy) + (tMax - deltaCutInEnergy) * ( (5. + 1. /var)*0.25 / (energy*energy) -
								    beta2 / (tMax * deltaCutInEnergy) ) / 3. ;
	}
      totalCrossSection *= twopi_mc2_rcl2 * atomicNumber / beta2 ;
    }

  //std::cout << "Microscopic = " << totalCrossSection/barn 
  //    << ", e = " << kineticEnergy/MeV <<std:: endl; 

  return totalCrossSection ;
}



// --------------------------------------------------------------------
G4double G4hImpactIonisation::GetMeanFreePath(const G4Track& track,
					      G4double, // previousStepSize
					      enum G4ForceCondition* condition)
{
  const G4DynamicParticle* dynamicParticle = track.GetDynamicParticle();
  const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
  const G4Material* material = couple->GetMaterial();

  G4double meanFreePath = DBL_MAX;
  // ---- MGP ---- What is the meaning of the local variable isOutOfRange?
  G4bool isOutRange = false;

  *condition = NotForced ;

  G4double kineticEnergy = (dynamicParticle->GetKineticEnergy())*initialMass/(dynamicParticle->GetMass());
  charge = dynamicParticle->GetCharge()/eplus;
  chargeSquare = theIonEffChargeModel->TheValue(dynamicParticle, material);

  if (kineticEnergy < LowestKineticEnergy) 
    {
      meanFreePath = DBL_MAX;
    }
  else 
    {
      if (kineticEnergy > HighestKineticEnergy)  kineticEnergy = HighestKineticEnergy;
      meanFreePath = (((*theMeanFreePathTable)(couple->GetIndex()))->
		      GetValue(kineticEnergy,isOutRange))/chargeSquare;
    }

  return meanFreePath ;
}


// --------------------------------------------------------------------
G4double G4hImpactIonisation::GetConstraints(const G4DynamicParticle* particle,
					     const G4MaterialCutsCouple* couple)
{
  // returns the Step limit
  // dEdx is calculated as well as the range
  // based on Effective Charge Approach

  const G4Material* material = couple->GetMaterial();
  G4Proton* proton = G4Proton::Proton();
  G4AntiProton* antiproton = G4AntiProton::AntiProton();

  G4double stepLimit = 0.;
  G4double dx, highEnergy;

  G4double massRatio = proton_mass_c2/(particle->GetMass()) ;
  G4double kineticEnergy = particle->GetKineticEnergy() ;

  // Scale the kinetic energy

  G4double tScaled = kineticEnergy*massRatio ;
  fBarkas = 0.;

  if (charge > 0.) 
    {
      highEnergy = protonHighEnergy ;
      fRangeNow = G4EnergyLossTables::GetRange(proton, tScaled, couple);
      dx = G4EnergyLossTables::GetRange(proton, highEnergy, couple);
      fdEdx = G4EnergyLossTables::GetDEDX(proton, tScaled, couple)
	* chargeSquare ;
      
      // Correction for positive ions
      if (theBarkas && tScaled > highEnergy) 
	{ 
	  fBarkas = BarkasTerm(material,tScaled)*std::sqrt(chargeSquare)*chargeSquare
	    + BlochTerm(material,tScaled,chargeSquare);
	}
      // Antiprotons and negative hadrons
    } 
  else 
    {
      highEnergy = antiprotonHighEnergy ;
      fRangeNow = G4EnergyLossTables::GetRange(antiproton, tScaled, couple);
      dx = G4EnergyLossTables::GetRange(antiproton, highEnergy, couple);
      fdEdx = G4EnergyLossTables::GetDEDX(antiproton, tScaled, couple) * chargeSquare ;
      
      if (theBarkas && tScaled > highEnergy) 
	{ 
	  fBarkas = -BarkasTerm(material,tScaled)*std::sqrt(chargeSquare)*chargeSquare
	    + BlochTerm(material,tScaled,chargeSquare);
	} 
    } 
  /*
    const G4Material* mat = couple->GetMaterial();
    G4double fac = gram/(MeV*cm2*mat->GetDensity());
    G4cout << particle->GetDefinition()->GetParticleName()
    << " in " << mat->GetName()
    << " E(MeV)= " << kineticEnergy/MeV
    << " dedx(MeV*cm^2/g)= " << fdEdx*fac
    << " barcas(MeV*cm^2/gram)= " << fBarkas*fac
    << " Q^2= " << chargeSquare
    << G4endl;
  */
  // scaling back
  fRangeNow /= (chargeSquare*massRatio) ;
  dx        /= (chargeSquare*massRatio) ;

  stepLimit  = fRangeNow ;
  G4double r = std::min(finalRange, couple->GetProductionCuts()
			->GetProductionCut(idxG4ElectronCut));

  if (fRangeNow > r) 
    {  
      stepLimit = dRoverRange*fRangeNow + r*(1.0 - dRoverRange)*(2.0 - r/fRangeNow);
      if (stepLimit > fRangeNow) stepLimit = fRangeNow;
    }  
  //   compute the (random) Step limit in standard energy range
  if(tScaled > highEnergy ) 
    {    
      // add Barkas correction directly to dedx
      fdEdx  += fBarkas;
      
      if(stepLimit > fRangeNow - dx*0.9) stepLimit = fRangeNow - dx*0.9 ;
      
      // Step limit in low energy range
    } 
  else 
    {
      G4double x = dx*paramStepLimit;
      if (stepLimit > x) stepLimit = x;
    }
  return stepLimit;
}


// --------------------------------------------------------------------
G4VParticleChange* G4hImpactIonisation::AlongStepDoIt(const G4Track& track,
						      const G4Step& step)
{
  // compute the energy loss after a step
  G4Proton* proton = G4Proton::Proton();
  G4AntiProton* antiproton = G4AntiProton::AntiProton();
  G4double finalT = 0.;

  aParticleChange.Initialize(track) ;

  const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
  const G4Material* material = couple->GetMaterial();

  // get the actual (true) Step length from step
  const G4double stepLength = step.GetStepLength() ;

  const G4DynamicParticle* particle = track.GetDynamicParticle() ;

  G4double kineticEnergy = particle->GetKineticEnergy() ;
  G4double massRatio = proton_mass_c2/(particle->GetMass()) ;
  G4double tScaled = kineticEnergy * massRatio ;
  G4double eLoss = 0.0 ;
  G4double nLoss = 0.0 ;


  // very small particle energy
  if (kineticEnergy < MinKineticEnergy) 
    {
      eLoss = kineticEnergy ;
      // particle energy outside tabulated energy range
    } 
  
  else if( kineticEnergy > HighestKineticEnergy) 
    {
      eLoss = stepLength * fdEdx ;
      // big step
    } 
  else if (stepLength >= fRangeNow ) 
    {
      eLoss = kineticEnergy ;
      
      // tabulated range
    } 
  else 
    {
      // step longer than linear step limit
      if(stepLength > linLossLimit * fRangeNow) 
	{
	  G4double rScaled = fRangeNow * massRatio * chargeSquare ;
	  G4double sScaled = stepLength * massRatio * chargeSquare ;
	  
	  if(charge > 0.0) 
	    {
	      eLoss = G4EnergyLossTables::GetPreciseEnergyFromRange(proton,rScaled, couple) -
		G4EnergyLossTables::GetPreciseEnergyFromRange(proton,rScaled-sScaled,couple) ;
	    
	    } 
	  else 
	    {
	      // Antiproton
	      eLoss = G4EnergyLossTables::GetPreciseEnergyFromRange(antiproton,rScaled,couple) -
		G4EnergyLossTables::GetPreciseEnergyFromRange(antiproton,rScaled-sScaled,couple) ;
	    }
	  eLoss /= massRatio ;
	  
	  // Barkas correction at big step      
	  eLoss += fBarkas * stepLength;
	  
	  // step shorter than linear step limit
	} 
      else 
	{
	  eLoss = stepLength *fdEdx  ;
	}
      if (nStopping && tScaled < protonHighEnergy) 
	{
	  nLoss = (theNuclearStoppingModel->TheValue(particle, material)) * stepLength;
	}
    }
  
  if (eLoss < 0.0) eLoss = 0.0;

  finalT = kineticEnergy - eLoss - nLoss;

  if ( EnlossFlucFlag && 0.0 < eLoss && finalT > MinKineticEnergy) 
    {
      
      //  now the electron loss with fluctuation
      eLoss = ElectronicLossFluctuation(particle, couple, eLoss, stepLength) ;
      if (eLoss < 0.0) eLoss = 0.0;
      finalT = kineticEnergy - eLoss - nLoss;
    }

  //  stop particle if the kinetic energy <= MinKineticEnergy
  if (finalT*massRatio <= MinKineticEnergy ) 
    {
      
      finalT = 0.0;
      if (!particle->GetDefinition()->GetProcessManager()->GetAtRestProcessVector()->size())
        aParticleChange.ProposeTrackStatus(fStopAndKill);
      else
        aParticleChange.ProposeTrackStatus(fStopButAlive);
    }

  aParticleChange.ProposeEnergy( finalT );
  eLoss = kineticEnergy-finalT;

  aParticleChange.ProposeLocalEnergyDeposit(eLoss);
  return &aParticleChange ;
}



// --------------------------------------------------------------------
G4double G4hImpactIonisation::ProtonParametrisedDEDX(const G4MaterialCutsCouple* couple,
						     G4double kineticEnergy) const
{
  const G4Material* material = couple->GetMaterial();
  G4Proton* proton = G4Proton::Proton();
  G4double eLoss = 0.;

  // Free Electron Gas Model
  if(kineticEnergy < protonLowEnergy) {
    eLoss = (protonModel->TheValue(proton, material, protonLowEnergy))
      * std::sqrt(kineticEnergy/protonLowEnergy) ;

    // Parametrisation
  } else {
    eLoss = protonModel->TheValue(proton, material, kineticEnergy) ;
  }

  // Delta rays energy
  eLoss -= DeltaRaysEnergy(couple,kineticEnergy,proton_mass_c2) ;

  if(verboseLevel > 2) {
    G4cout << "p E(MeV)= " << kineticEnergy/MeV
           << " dE/dx(MeV/mm)= " << eLoss*mm/MeV
           << " for " << material->GetName()
           << " model: " << protonModel << G4endl;
  }

  if(eLoss < 0.0) eLoss = 0.0 ;

  return eLoss ;
}



// --------------------------------------------------------------------
G4double G4hImpactIonisation::AntiProtonParametrisedDEDX(const G4MaterialCutsCouple* couple,
							 G4double kineticEnergy) const
{
  const G4Material* material = couple->GetMaterial();
  G4AntiProton* antiproton = G4AntiProton::AntiProton();
  G4double eLoss = 0.0 ;

  // Antiproton model is used
  if(antiprotonModel->IsInCharge(antiproton,material)) {
    if(kineticEnergy < antiprotonLowEnergy) {
      eLoss = antiprotonModel->TheValue(antiproton,material,antiprotonLowEnergy)
	* std::sqrt(kineticEnergy/antiprotonLowEnergy) ;

      // Parametrisation
    } else {
      eLoss = antiprotonModel->TheValue(antiproton,material,
					kineticEnergy);
    }

    // The proton model is used + Barkas correction
  } else {
    if(kineticEnergy < protonLowEnergy) {
      eLoss = protonModel->TheValue(G4Proton::Proton(),material,protonLowEnergy)
	* std::sqrt(kineticEnergy/protonLowEnergy) ;

      // Parametrisation
    } else {
      eLoss = protonModel->TheValue(G4Proton::Proton(),material,
				    kineticEnergy);
    }
    //if(theBarkas) eLoss -= 2.0*BarkasTerm(material, kineticEnergy);
  }

  // Delta rays energy
  eLoss -= DeltaRaysEnergy(couple,kineticEnergy,proton_mass_c2) ;

  if(verboseLevel > 2) {
    G4cout << "pbar E(MeV)= " << kineticEnergy/MeV
           << " dE/dx(MeV/mm)= " << eLoss*mm/MeV
           << " for " << material->GetName()
           << " model: " << protonModel << G4endl;
  }

  if(eLoss < 0.0) eLoss = 0.0 ;

  return eLoss ;
}


// --------------------------------------------------------------------
G4double G4hImpactIonisation::DeltaRaysEnergy(const G4MaterialCutsCouple* couple,
					      G4double kineticEnergy,
					      G4double particleMass) const
{
  G4double dLoss = 0.;

  G4double deltaCutNow = cutForDelta[(couple->GetIndex())] ;
  const G4Material* material = couple->GetMaterial();
  G4double electronDensity = material->GetElectronDensity();
  G4double excitationEnergy = material->GetIonisation()->GetMeanExcitationEnergy();

  G4double tau = kineticEnergy / particleMass ;
  G4double rateMass = electron_mass_c2/particleMass ;

  // some local variables

  G4double gamma = tau + 1.0 ;
  G4double bg2 = tau*(tau+2.0) ;
  G4double beta2 = bg2/(gamma*gamma) ;
  G4double tMax = 2.*electron_mass_c2*bg2/(1.0+2.0*gamma*rateMass+rateMass*rateMass) ;

  // Validity range for delta electron cross section
  G4double deltaCut = std::max(deltaCutNow, excitationEnergy);

  if ( deltaCut < tMax) 
    {
      G4double x = deltaCut / tMax ;
      dLoss = ( beta2 * (x-1.) - std::log(x) ) * twopi_mc2_rcl2 * electronDensity / beta2 ;
    }
  return dLoss ;
}


// -------------------------------------------------------------------------

G4VParticleChange* G4hImpactIonisation::PostStepDoIt(const G4Track& track,
						     const G4Step& step)
{
  // Units are expressed in GEANT4 internal units.

  //  std::cout << "----- Calling PostStepDoIt ----- " << std::endl;

  aParticleChange.Initialize(track) ;
  const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();  
  const G4DynamicParticle* aParticle = track.GetDynamicParticle() ;
  
  // Some kinematics

  G4ParticleDefinition* definition = track.GetDefinition();
  G4double mass = definition->GetPDGMass();
  G4double kineticEnergy = aParticle->GetKineticEnergy();
  G4double totalEnergy = kineticEnergy + mass ;
  G4double pSquare = kineticEnergy *( totalEnergy + mass) ;
  G4double eSquare = totalEnergy * totalEnergy;
  G4double betaSquare = pSquare / eSquare;
  G4ThreeVector particleDirection = aParticle->GetMomentumDirection() ;

  G4double gamma = kineticEnergy / mass + 1.;
  G4double r = electron_mass_c2 / mass;
  G4double tMax = 2. * electron_mass_c2 *(gamma*gamma - 1.) / (1. + 2.*gamma*r + r*r);

  // Validity range for delta electron cross section
  G4double deltaCut = cutForDelta[couple->GetIndex()];

  // This should not be a case
  if (deltaCut >= tMax)
    return G4VContinuousDiscreteProcess::PostStepDoIt(track,step);

  G4double xc = deltaCut / tMax;
  G4double rate = tMax / totalEnergy;
  rate = rate*rate ;
  G4double spin = aParticle->GetDefinition()->GetPDGSpin() ;

  // Sampling follows ...
  G4double x = 0.;
  G4double gRej = 0.;

  do {
    x = xc / (1. - (1. - xc) * G4UniformRand());
    
    if (0.0 == spin) 
      {
	gRej = 1.0 - betaSquare * x ;
      }
    else if (0.5 == spin) 
      {
	gRej = (1. - betaSquare * x + 0.5 * x*x * rate) / (1. + 0.5 * rate) ;
      } 
    else 
      {
	gRej = (1. - betaSquare * x ) * (1. + x/(3.*xc)) +
	  x*x * rate * (1. + 0.5*x/xc) / 3.0 /
	  (1. + 1./(3.*xc) + rate *(1.+ 0.5/xc) / 3.);
      }
    
  } while( G4UniformRand() > gRej );

  G4double deltaKineticEnergy = x * tMax;
  G4double deltaTotalMomentum = std::sqrt(deltaKineticEnergy * 
					  (deltaKineticEnergy + 2. * electron_mass_c2 ));
  G4double totalMomentum = std::sqrt(pSquare) ;
  G4double cosTheta = deltaKineticEnergy * (totalEnergy + electron_mass_c2) / (deltaTotalMomentum*totalMomentum);

  //  protection against cosTheta > 1 or < -1 
  if ( cosTheta < -1. ) cosTheta = -1.;
  if ( cosTheta > 1. ) cosTheta = 1.;

  //  direction of the delta electron 
  G4double phi = twopi * G4UniformRand();
  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  G4double dirX = sinTheta * std::cos(phi);
  G4double dirY = sinTheta * std::sin(phi);
  G4double dirZ = cosTheta;

  G4ThreeVector deltaDirection(dirX,dirY,dirZ);
  deltaDirection.rotateUz(particleDirection);

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* deltaRay = new G4DynamicParticle;
  deltaRay->SetKineticEnergy( deltaKineticEnergy );
  deltaRay->SetMomentumDirection(deltaDirection.x(),
				 deltaDirection.y(),
				 deltaDirection.z());
  deltaRay->SetDefinition(G4Electron::Electron());

  // fill aParticleChange
  G4double finalKineticEnergy = kineticEnergy - deltaKineticEnergy;
  std::size_t totalNumber = 1;

  // Atomic relaxation

  // ---- MGP ---- Temporary limitation: currently PIXE only for incident protons

  std::size_t nSecondaries = 0;
  std::vector<G4DynamicParticle*>* secondaryVector = 0;

  if (definition == G4Proton::ProtonDefinition())
    {
      const G4Material* material = couple->GetMaterial();

      // Lazy initialization of pixeCrossSectionHandler
      if (pixeCrossSectionHandler == 0)
	{
	  // Instantiate pixeCrossSectionHandler with selected shell cross section models
	  // Ownership of interpolation is transferred to pixeCrossSectionHandler
	  G4IInterpolator* interpolation = new G4LogLogInterpolator();
	  pixeCrossSectionHandler = 
	    new G4PixeCrossSectionHandler(interpolation,modelK,modelL,modelM,eMinPixe,eMaxPixe);
	  G4String fileName("proton");
	  pixeCrossSectionHandler->LoadShellData(fileName);
	  //	  pixeCrossSectionHandler->PrintData();
	}
      
      // Select an atom in the current material based on the total shell cross sections
      G4int Z = pixeCrossSectionHandler->SelectRandomAtom(material,kineticEnergy);
      //      std::cout << "G4hImpactIonisation::PostStepDoIt - Z = " << Z << std::endl;

      //      G4double microscopicCross = MicroscopicCrossSection(*definition,
      //       			         	  kineticEnergy,
      //					  Z, deltaCut);
  //    G4double crossFromShells = pixeCrossSectionHandler->FindValue(Z,kineticEnergy);

      //std::cout << "G4hImpactIonisation: Z= "
      //		<< Z
      //		<< ", energy = "
      //		<< kineticEnergy/MeV
      //		<<" MeV, microscopic = "
      //		<< microscopicCross/barn 
      //		<< " barn, from shells = "
      //		<< crossFromShells/barn
      //		<< " barn" 
      //		<< std::endl;

      // Select a shell in the target atom based on the individual shell cross sections
      G4int shellIndex = pixeCrossSectionHandler->SelectRandomShell(Z,kineticEnergy);
    
      G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
      const G4AtomicShell* atomicShell = transitionManager->Shell(Z,shellIndex);
      G4double bindingEnergy = atomicShell->BindingEnergy();
      
      //     if (verboseLevel > 1) 
      //       {
      //	 G4cout << "G4hImpactIonisation::PostStepDoIt - Z = " 
      //         << Z 
      //	 << ", shell = " 
      //	 << shellIndex
      //	 << ", bindingE (keV) = " 
      //	 << bindingEnergy/keV
      //	 << G4endl;
      //       }
      
      // Generate PIXE if binding energy larger than cut for photons or electrons
      
      G4ParticleDefinition* type = 0;
      
      if (finalKineticEnergy >= bindingEnergy)
	//	  && (bindingEnergy >= minGammaEnergy ||  bindingEnergy >= minElectronEnergy) ) 
	{
	  // Vacancy in subshell shellIndex; shellId is the subshell identifier in EADL jargon 
	  G4int shellId = atomicShell->ShellId();
	  // Atomic relaxation: generate secondaries
	  secondaryVector = atomicDeexcitation.GenerateParticles(Z, shellId);

	  // ---- Debug ----
	  //std::cout << "ShellId = "
	  //    <<shellId << " ---- Atomic relaxation secondaries: ---- " 
	  //    << secondaryVector->size()
	  //    << std::endl;

	  // ---- End debug ---

	  if (secondaryVector != 0) 
	    {
	      nSecondaries = secondaryVector->size();
	      for (std::size_t i = 0; i<nSecondaries; i++) 
		{
		  G4DynamicParticle* aSecondary = (*secondaryVector)[i];
		  if (aSecondary) 
		    {
		      G4double e = aSecondary->GetKineticEnergy();
		      type = aSecondary->GetDefinition();

		      // ---- Debug ----
		      //if (type == G4Gamma::GammaDefinition())
		      //	{			
		      //	  std::cout << "Z = " << Z 
		      //		    << ", shell: " << shellId
		      //		    << ", PIXE photon energy (keV) = " << e/keV 
		      //		    << std::endl;
		      //	}
		      // ---- End debug ---

		      if (e < finalKineticEnergy &&
			  ((type == G4Gamma::Gamma() && e > minGammaEnergy ) ||
			   (type == G4Electron::Electron() && e > minElectronEnergy ))) 
			{
			  // Subtract the energy of the emitted secondary from the primary
			  finalKineticEnergy -= e;
			  totalNumber++;
			  // ---- Debug ----
			  //if (type == G4Gamma::GammaDefinition())
			  //	{			
			  //	  std::cout << "Z = " << Z 
			  //		    << ", shell: " << shellId
			  //		    << ", PIXE photon energy (keV) = " << e/keV
			  //		    << std::endl;
			  //	}
			  // ---- End debug ---
			} 
		      else 
			{
			  // The atomic relaxation product has energy below the cut
			  // ---- Debug ----
			  // if (type == G4Gamma::GammaDefinition())
			  //	{			
			  //	  std::cout << "Z = " << Z 
			  //
			  //		    << ", PIXE photon energy = " << e/keV 
			  //  		    << " keV below threshold " << minGammaEnergy/keV << " keV"
			  //		    << std::endl;
			  //	}   
			  // ---- End debug ---
     
			  delete aSecondary;
			  (*secondaryVector)[i] = 0;
			}
		    }
		}
	    }
	}
    }


  // Save delta-electrons

  G4double eDeposit = 0.;

  if (finalKineticEnergy > MinKineticEnergy)
    {
      G4double finalPx = totalMomentum*particleDirection.x() - deltaTotalMomentum*deltaDirection.x();
      G4double finalPy = totalMomentum*particleDirection.y() - deltaTotalMomentum*deltaDirection.y();
      G4double finalPz = totalMomentum*particleDirection.z() - deltaTotalMomentum*deltaDirection.z();
      G4double finalMomentum = std::sqrt(finalPx*finalPx + finalPy*finalPy + finalPz*finalPz) ;
      finalPx /= finalMomentum;
      finalPy /= finalMomentum;
      finalPz /= finalMomentum;

      aParticleChange.ProposeMomentumDirection( finalPx,finalPy,finalPz );
    }
  else
    {
      eDeposit = finalKineticEnergy;
      finalKineticEnergy = 0.;
      aParticleChange.ProposeMomentumDirection(particleDirection.x(),
					       particleDirection.y(),
					       particleDirection.z());
      if(!aParticle->GetDefinition()->GetProcessManager()->
	 GetAtRestProcessVector()->size())
        aParticleChange.ProposeTrackStatus(fStopAndKill);
      else
        aParticleChange.ProposeTrackStatus(fStopButAlive);
    }

  aParticleChange.ProposeEnergy(finalKineticEnergy);
  aParticleChange.ProposeLocalEnergyDeposit (eDeposit);
  aParticleChange.SetNumberOfSecondaries((G4int)totalNumber);
  aParticleChange.AddSecondary(deltaRay);

  // ---- Debug ----
  //  std::cout << "RDHadronIonisation - finalKineticEnergy (MeV) = " 
  //	    << finalKineticEnergy/MeV 
  //	    << ", delta KineticEnergy (keV) = " 
  //	    << deltaKineticEnergy/keV 
  //	    << ", energy deposit (MeV) = "
  //	    << eDeposit/MeV
  //	    << std::endl;
  // ---- End debug ---
  
  // Save Fluorescence and Auger

  if (secondaryVector != 0) 
    { 
      for (std::size_t l = 0; l < nSecondaries; l++) 
	{ 
	  G4DynamicParticle* secondary = (*secondaryVector)[l];
	  if (secondary) aParticleChange.AddSecondary(secondary);

	  // ---- Debug ----
	  //if (secondary != 0) 
	  // {
	  //   if (secondary->GetDefinition() == G4Gamma::GammaDefinition())
	  //	{
	  //	  G4double eX = secondary->GetKineticEnergy();			
	  //	  std::cout << " PIXE photon of energy " << eX/keV 
	  //		    << " keV added to ParticleChange; total number of secondaries is " << totalNumber 
	  //		    << std::endl;
	  //	}
	  //}
	  // ---- End debug ---	  

	}
      delete secondaryVector;
    }

  return G4VContinuousDiscreteProcess::PostStepDoIt(track,step);
}

// -------------------------------------------------------------------------

G4double G4hImpactIonisation::ComputeDEDX(const G4ParticleDefinition* aParticle,
					  const G4MaterialCutsCouple* couple,
					  G4double kineticEnergy)
{
  const G4Material* material = couple->GetMaterial();
  G4Proton* proton = G4Proton::Proton();
  G4AntiProton* antiproton = G4AntiProton::AntiProton();
  G4double dedx = 0.;

  G4double tScaled = kineticEnergy * proton_mass_c2 / (aParticle->GetPDGMass()) ;
  charge  = aParticle->GetPDGCharge() ;

  if( charge > 0.) 
    {
      if (tScaled > protonHighEnergy) 
	{
	  dedx = G4EnergyLossTables::GetDEDX(proton,tScaled,couple) ;  
	}
      else 
	{
	  dedx = ProtonParametrisedDEDX(couple,tScaled) ;
	}
    } 
  else 
    {
      if (tScaled > antiprotonHighEnergy) 
	{
	  dedx = G4EnergyLossTables::GetDEDX(antiproton,tScaled,couple);
	} 
      else
	{
	  dedx = AntiProtonParametrisedDEDX(couple,tScaled) ;
	}
    }
  dedx *= theIonEffChargeModel->TheValue(aParticle, material, kineticEnergy) ;
  
  return dedx ;
}


G4double G4hImpactIonisation::BarkasTerm(const G4Material* material,
					 G4double kineticEnergy) const
//Function to compute the Barkas term for protons:
//
//Ref. Z_1^3 effect in the stopping power of matter for charged particles
//     J.C Ashley and R.H.Ritchie
//     Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397
//
{
  static G4ThreadLocal G4double FTable[47][2] = {
    { 0.02, 21.5},
    { 0.03, 20.0},
    { 0.04, 18.0},
    { 0.05, 15.6},
    { 0.06, 15.0},
    { 0.07, 14.0},
    { 0.08, 13.5},
    { 0.09, 13.},
    { 0.1,  12.2},
    { 0.2,  9.25},
    { 0.3,  7.0},
    { 0.4,  6.0},
    { 0.5,  4.5},
    { 0.6,  3.5},
    { 0.7,  3.0},
    { 0.8,  2.5},
    { 0.9,  2.0},
    { 1.0,  1.7},
    { 1.2,  1.2},
    { 1.3,  1.0},
    { 1.4,  0.86},
    { 1.5,  0.7},
    { 1.6,  0.61},
    { 1.7,  0.52},
    { 1.8,  0.5},
    { 1.9,  0.43},
    { 2.0,  0.42},
    { 2.1,  0.3},
    { 2.4,  0.2},
    { 3.0,  0.13},
    { 3.08, 0.1},
    { 3.1,  0.09},
    { 3.3,  0.08},
    { 3.5,  0.07},
    { 3.8,  0.06},
    { 4.0,  0.051},
    { 4.1,  0.04},
    { 4.8,  0.03},
    { 5.0,  0.024},
    { 5.1,  0.02},
    { 6.0,  0.013},
    { 6.5,  0.01},
    { 7.0,  0.009},
    { 7.1,  0.008},
    { 8.0,  0.006},
    { 9.0,  0.0032},
    { 10.0, 0.0025} };

  // Information on particle and material
  G4double kinE  = kineticEnergy ;
  if(0.5*MeV > kinE) kinE = 0.5*MeV ;
  G4double gamma = 1.0 + kinE / proton_mass_c2 ;
  G4double beta2 = 1.0 - 1.0/(gamma*gamma) ;
  if(0.0 >= beta2) return 0.0;

  G4double BTerm = 0.0;
  //G4double AMaterial = 0.0;
  G4double ZMaterial = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int numberOfElements = (G4int)material->GetNumberOfElements();

  for (G4int i = 0; i<numberOfElements; i++) {

    //AMaterial = (*theElementVector)[i]->GetA()*mole/g;
    ZMaterial = (*theElementVector)[i]->GetZ();

    G4double X = 137.0 * 137.0 * beta2 / ZMaterial;

    // Variables to compute L_1
    G4double Eta0Chi = 0.8;
    G4double EtaChi = Eta0Chi * ( 1.0 + 6.02*std::pow( ZMaterial,-1.19 ) );
    G4double W = ( EtaChi * std::pow( ZMaterial,1.0/6.0 ) ) / std::sqrt(X);
    G4double FunctionOfW = FTable[46][1]*FTable[46][0]/W ;

    for(G4int j=0; j<47; j++) {

      if( W < FTable[j][0] ) {

        if(0 == j) {
     	  FunctionOfW = FTable[0][1] ;

        } else {
          FunctionOfW = (FTable[j][1] - FTable[j-1][1]) * (W - FTable[j-1][0])
	    / (FTable[j][0] - FTable[j-1][0])
	    +  FTable[j-1][1] ;
	}

        break;
      }

    }

    BTerm += FunctionOfW /( std::sqrt(ZMaterial * X) * X);
  }

  BTerm *= twopi_mc2_rcl2 * (material->GetElectronDensity()) / beta2 ;

  return BTerm;
}



G4double G4hImpactIonisation::BlochTerm(const G4Material* material,
					G4double kineticEnergy,
					G4double cSquare) const
//Function to compute the Bloch term for protons:
//
//Ref. Z_1^3 effect in the stopping power of matter for charged particles
//     J.C Ashley and R.H.Ritchie
//     Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397
//
{
  G4double eLoss = 0.0 ;
  G4double gamma = 1.0 + kineticEnergy / proton_mass_c2 ;
  G4double beta2 = 1.0 - 1.0/(gamma*gamma) ;
  G4double y = cSquare / (137.0*137.0*beta2) ;

  if(y < 0.05) {
    eLoss = 1.202 ;

  } else {
    eLoss = 1.0 / (1.0 + y) ;
    G4double de = eLoss ;

    for(G4int i=2; de>eLoss*0.01; i++) {
      de = 1.0/( i * (i*i + y)) ;
      eLoss += de ;
    }
  }
  eLoss *= -1.0 * y * cSquare * twopi_mc2_rcl2 *
    (material->GetElectronDensity()) / beta2 ;

  return eLoss;
}



G4double G4hImpactIonisation::ElectronicLossFluctuation(
							const G4DynamicParticle* particle,
							const G4MaterialCutsCouple* couple,
							G4double meanLoss,
							G4double step) const
//  calculate actual loss from the mean loss
//  The model used to get the fluctuation is essentially the same
// as in Glandz in Geant3.
{
  // data members to speed up the fluctuation calculation
  //  G4int imat ;
  //  G4double f1Fluct,f2Fluct,e1Fluct,e2Fluct,rateFluct,ipotFluct;
  //  G4double e1LogFluct,e2LogFluct,ipotLogFluct;

  static const G4double minLoss = 1.*eV ;
  static const G4double kappa = 10. ;
  static const G4double theBohrBeta2 = 50.0 * keV/proton_mass_c2 ;

  const G4Material* material = couple->GetMaterial();
  G4int    imaterial   = couple->GetIndex() ;
  G4double ipotFluct   = material->GetIonisation()->GetMeanExcitationEnergy() ;
  G4double electronDensity = material->GetElectronDensity() ;
  G4double zeff = electronDensity/(material->GetTotNbOfAtomsPerVolume()) ;

  // get particle data
  G4double tkin   = particle->GetKineticEnergy();
  G4double particleMass = particle->GetMass() ;
  G4double deltaCutInKineticEnergyNow = cutForDelta[imaterial];

  // shortcut for very very small loss
  if(meanLoss < minLoss) return meanLoss ;

  // Validity range for delta electron cross section
  G4double threshold = std::max(deltaCutInKineticEnergyNow,ipotFluct);
  G4double loss, siga;

  G4double rmass = electron_mass_c2/particleMass;
  G4double tau   = tkin/particleMass;
  G4double tau1 = tau+1.0;
  G4double tau2 = tau*(tau+2.);
  G4double tMax    = 2.*electron_mass_c2*tau2/(1.+2.*tau1*rmass+rmass*rmass);


  if(tMax > threshold) tMax = threshold;
  G4double beta2 = tau2/(tau1*tau1);

  // Gaussian fluctuation
  if(meanLoss > kappa*tMax || tMax < kappa*ipotFluct )
    {
      siga = tMax * (1.0-0.5*beta2) * step * twopi_mc2_rcl2
	* electronDensity / beta2 ;

      // High velocity or negatively charged particle
      if( beta2 > 3.0*theBohrBeta2*zeff || charge < 0.0) {
	siga = std::sqrt( siga * chargeSquare ) ;

	// Low velocity - additional ion charge fluctuations according to
	// Q.Yang et al., NIM B61(1991)149-155.
      } else {
	G4double chu = theIonChuFluctuationModel->TheValue(particle, material);
	G4double yang = theIonYangFluctuationModel->TheValue(particle, material);
	siga = std::sqrt( siga * (chargeSquare * chu + yang)) ;
      }

      do {
        loss = G4RandGauss::shoot(meanLoss,siga);
      } while (loss < 0. || loss > 2.0*meanLoss);
      return loss;
    }

  // Non Gaussian fluctuation
  static const G4double probLim = 0.01 ;
  static const G4double sumaLim = -std::log(probLim) ;
  static const G4double alim = 10.;

  G4double suma,w1,w2,C,e0,lossc,w;
  G4double a1,a2,a3;
  G4int p1,p2,p3;
  G4int nb;
  G4double corrfac, na,alfa,rfac,namean,sa,alfa1,ea,sea;
  G4double dp3;

  G4double f1Fluct     = material->GetIonisation()->GetF1fluct();
  G4double f2Fluct     = material->GetIonisation()->GetF2fluct();
  G4double e1Fluct     = material->GetIonisation()->GetEnergy1fluct();
  G4double e2Fluct     = material->GetIonisation()->GetEnergy2fluct();
  G4double e1LogFluct  = material->GetIonisation()->GetLogEnergy1fluct();
  G4double e2LogFluct  = material->GetIonisation()->GetLogEnergy2fluct();
  G4double rateFluct   = material->GetIonisation()->GetRateionexcfluct();
  G4double ipotLogFluct= material->GetIonisation()->GetLogMeanExcEnergy();

  w1 = tMax/ipotFluct;
  w2 = std::log(2.*electron_mass_c2*tau2);

  C = meanLoss*(1.-rateFluct)/(w2-ipotLogFluct-beta2);

  a1 = C*f1Fluct*(w2-e1LogFluct-beta2)/e1Fluct;
  a2 = C*f2Fluct*(w2-e2LogFluct-beta2)/e2Fluct;
  a3 = rateFluct*meanLoss*(tMax-ipotFluct)/(ipotFluct*tMax*std::log(w1));
  if(a1 < 0.0) a1 = 0.0;
  if(a2 < 0.0) a2 = 0.0;
  if(a3 < 0.0) a3 = 0.0;

  suma = a1+a2+a3;

  loss = 0.;


  if(suma < sumaLim)             // very small Step
    {
      e0 = material->GetIonisation()->GetEnergy0fluct();

      if(tMax == ipotFluct)
	{
	  a3 = meanLoss/e0;

	  if(a3>alim)
	    {
	      siga=std::sqrt(a3) ;
	      p3 = std::max(0,G4int(G4RandGauss::shoot(a3,siga)+0.5));
	    }
	  else
	    p3 = (G4int)G4Poisson(a3);

	  loss = p3*e0 ;

	  if(p3 > 0)
	    loss += (1.-2.*G4UniformRand())*e0 ;

	}
      else
	{
	  tMax = tMax-ipotFluct+e0 ;
	  a3 = meanLoss*(tMax-e0)/(tMax*e0*std::log(tMax/e0));

	  if(a3>alim)
	    {
	      siga=std::sqrt(a3) ;
	      p3 = std::max(0,int(G4RandGauss::shoot(a3,siga)+0.5));
	    }
	  else
	    p3 = (G4int)G4Poisson(a3);

	  if(p3 > 0)
	    {
	      w = (tMax-e0)/tMax ;
	      if(p3 > nmaxCont2)
		{
		  dp3 = G4float(p3) ;
		  corrfac = dp3/G4float(nmaxCont2) ;
		  p3 = G4int(nmaxCont2) ;
		}
	      else
		corrfac = 1. ;

	      for(G4int i=0; i<p3; i++) loss += 1./(1.-w*G4UniformRand()) ;
	      loss *= e0*corrfac ;
	    }
	}
    }

  else                              // not so small Step
    {
      // excitation type 1
      if(a1>alim)
	{
	  siga=std::sqrt(a1) ;
	  p1 = std::max(0,G4int(G4RandGauss::shoot(a1,siga)+0.5));
	}
      else
	p1 = (G4int)G4Poisson(a1);

      // excitation type 2
      if(a2>alim)
	{
	  siga=std::sqrt(a2) ;
	  p2 = std::max(0,G4int(G4RandGauss::shoot(a2,siga)+0.5));
	}
      else
        p2 = (G4int)G4Poisson(a2);

      loss = p1*e1Fluct+p2*e2Fluct;

      // smearing to avoid unphysical peaks
      if(p2 > 0)
        loss += (1.-2.*G4UniformRand())*e2Fluct;
      else if (loss>0.)
        loss += (1.-2.*G4UniformRand())*e1Fluct;

      // ionisation .......................................
      if(a3 > 0.)
	{
	  if(a3>alim)
	    {
	      siga=std::sqrt(a3) ;
	      p3 = std::max(0,G4int(G4RandGauss::shoot(a3,siga)+0.5));
	    }
	  else
	    p3 = (G4int)G4Poisson(a3);

	  lossc = 0.;
	  if(p3 > 0)
	    {
	      na = 0.;
	      alfa = 1.;
	      if (p3 > nmaxCont2)
		{
		  dp3        = G4float(p3);
		  rfac       = dp3/(G4float(nmaxCont2)+dp3);
		  namean     = G4float(p3)*rfac;
		  sa         = G4float(nmaxCont1)*rfac;
		  na         = G4RandGauss::shoot(namean,sa);
		  if (na > 0.)
		    {
		      alfa   = w1*G4float(nmaxCont2+p3)/
			(w1*G4float(nmaxCont2)+G4float(p3));
		      alfa1  = alfa*std::log(alfa)/(alfa-1.);
		      ea     = na*ipotFluct*alfa1;
		      sea    = ipotFluct*std::sqrt(na*(alfa-alfa1*alfa1));
		      lossc += G4RandGauss::shoot(ea,sea);
		    }
		}

	      nb = G4int(G4float(p3)-na);
	      if (nb > 0)
		{
		  w2 = alfa*ipotFluct;
		  w  = (tMax-w2)/tMax;
		  for (G4int k=0; k<nb; k++) lossc += w2/(1.-w*G4UniformRand());
		}
	    }
	  loss += lossc;
	}
    }

  return loss ;
}



void G4hImpactIonisation::SetCutForSecondaryPhotons(G4double cut)
{
  minGammaEnergy = cut;
}



void G4hImpactIonisation::SetCutForAugerElectrons(G4double cut)
{
  minElectronEnergy = cut;
}



void G4hImpactIonisation::ActivateAugerElectronProduction(G4bool val)
{
  atomicDeexcitation.ActivateAugerElectronProduction(val);
}



void G4hImpactIonisation::PrintInfoDefinition() const
{
  G4String comments = "  Knock-on electron cross sections . ";
  comments += "\n        Good description above the mean excitation energy.\n";
  comments += "        Delta ray energy sampled from  differential Xsection.";

  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << LowestKineticEnergy / eV << " eV "
         << " to " << HighestKineticEnergy / TeV << " TeV "
         << " in " << TotBin << " bins."
	 << "\n        Electronic stopping power model is  "
	 << protonTable
         << "\n        from " << protonLowEnergy / keV << " keV "
         << " to " << protonHighEnergy / MeV << " MeV " << "." << G4endl ;
  G4cout << "\n        Parametrisation model for antiprotons is  "
         << antiprotonTable
         << "\n        from " << antiprotonLowEnergy / keV << " keV "
         << " to " << antiprotonHighEnergy / MeV << " MeV " << "." << G4endl ;
  if(theBarkas){
    G4cout << "        Parametrization of the Barkas effect is switched on."
	   << G4endl ;
  }
  if(nStopping) {
    G4cout << "        Nuclear stopping power model is " << theNuclearTable
	   << G4endl ;
  }

  G4bool printHead = true;

  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  // loop for materials

  for (G4int j=0 ; j < numOfCouples; ++j) {

    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);
    const G4Material* material= couple->GetMaterial();
    G4double deltaCutNow = cutForDelta[(couple->GetIndex())] ;
    G4double excitationEnergy = material->GetIonisation()->GetMeanExcitationEnergy();

    if(excitationEnergy > deltaCutNow) {
      if(printHead) {
        printHead = false ;

        G4cout << "           material       min.delta energy(keV) " << G4endl;
        G4cout << G4endl;
      }

      G4cout << std::setw(20) << material->GetName()
	     << std::setw(15) << excitationEnergy/keV << G4endl;
    }
  }
}
