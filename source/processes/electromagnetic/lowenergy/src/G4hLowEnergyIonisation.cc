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
//
// -------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4hLowEnergyIonisation physics process -------
//                by Vladimir Ivanchenko, 14 July 1999 
//                was made on the base of G4hIonisation class
//                developed by Laszlo Urban
// ************************************************************
// It is the extention of the ionisation process for the slow 
// charged hadrons.
// ************************************************************
// 28 July   1999 V.Ivanchenko cleen up
// 17 August 1999 G.Mancinelli added ICRU parametrisations for protons  
// 20 August 1999 G.Mancinelli added ICRU tables for alpha 
// 31 August 1999 V.Ivanchenko update and cleen up 
// 30 Sept.  1999 V.Ivanchenko minor upgrade 
// 12 Dec.   1999 S. Chauvie added Barkas correction 
// 19 Jan.   2000 V.Ivanchenko minor changing in Barkas corrections
// 02 April  2000 S. Chauvie linearization of Barkas effect
// 03 April  2000 V.Ivanchenko Nuclear Stopping power for antiprotons
// 23 May    2000 MG Pia  Clean up for QAO model 
// 24 May    2000 MG Pia  Code properly indented to improve legibility
// 17 July   2000 V.Ivanchenko Bug in scaling AlongStepDoIt method
// 25 July   2000 V.Ivanchenko New design iteration
// 17 August 2000 V.Ivanchenko Add ion fluctuation models
// 18 August 2000 V.Ivanchenko Bug fixed in GetConstrain
// 22 August 2000 V.Ivanchenko Insert paramStepLimit and
//                reorganise access to Barkas and Bloch terms  
// 04 Sept.  2000 V.Ivanchenko rename fluctuations
// 05 Sept.  2000 V.Ivanchenko clean up
// 03 Oct.   2000 V.Ivanchenko CodeWizard clean up
// 03 Nov.   2000 V.Ivanchenko MinKineticEnergy=LowestKineticEnergy=10eV
// 05 Nov.   2000 MG Pia - Removed const cast previously introduced to get
//                the code compiled (const G4Material* now introduced in 
//                electromagnetic/utils utils-V02-00-03 tag)
//                (this is going back and forth, to cope with Michel's
//                utils tag not being accepted yet by system testing)
// 21 Nov.  2000 V.Ivanchenko Fix a problem in fluctuations
// 23 Nov.  2000 V.Ivanchenko Ion type fluctuations only for charge>0
// 10 May   2001 V.Ivanchenko Clean up againist Linux compilation with -Wall
// 23 May   2001 V.Ivanchenko Minor fix in PostStepDoIt
// 07 June  2001 V.Ivanchenko Clean up AntiProtonDEDX + add print out
// 18 June  2001 V.Ivanchenko Cleanup print out
// 18 Oct.  2001 V.Ivanchenko Add fluorescence
// 30 Oct.  2001 V.Ivanchenko Add minGammaEnergy and minElectronEnergy
// 07 Dec   2001 V.Ivanchenko Add SetFluorescence method
// 15 Feb   2002 V.Ivanchenko Fix problem of Generic Ions
// 25 Mar   2002 V.Ivanchenko Fix problem of fluorescence below threshold  
// 28 Mar   2002 V.Ivanchenko Set fluorescence off by default
// 09 Apr   2002 V.Ivanchenko Fix table problem of GenericIons
// 28 May   2002 V.Ivanchenko Remove flag fStopAndKill
// 31 May   2002 V.Ivanchenko Add path of Fluo + Auger cuts to 
//                            AtomicDeexcitation
// 03 Jun   2002 MGP          Restore fStopAndKill
// 10 Jun   2002 V.Ivanchenko Restore fStopButAlive
// 12 Jun   2002 V.Ivanchenko Fix in fluctuations - if tmax<2*Ipot Gaussian
//                            fluctuations enables

// -----------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hLowEnergyIonisation.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4UnitsTable.hh"
#include "G4EnergyLossTables.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4AtomicDeexcitation.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4ShellVacancy.hh"
#include "G4hShellCrossSection.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4Gamma.hh"
#include "G4LogLogInterpolation.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4ProcessManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hLowEnergyIonisation::G4hLowEnergyIonisation(const G4String& processName)
  : G4hLowEnergyLoss(processName),
    theBetheBlochModel(0),
    theProtonModel(0),
    theAntiProtonModel(0),
    theIonEffChargeModel(0),
    theNuclearStoppingModel(0),
    theIonChuFluctuationModel(0),
    theIonYangFluctuationModel(0),
    theProtonTable("ICRU_R49p"),
    theAntiProtonTable("ICRU_R49p"),
    theNuclearTable("ICRU_R49"),
    nStopping(true),
    theBarkas(true),
    theMeanFreePathTable(0),
    paramStepLimit (0.005),
    shellVacancy(0),
    shellCS(0),
    theFluo(false)
{ 
  InitializeMe();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::InitializeMe()
{
  LowestKineticEnergy  = 10.0*eV ;
  HighestKineticEnergy = 100.0*TeV ;
  MinKineticEnergy     = 10.0*eV ; 
  TotBin               = 200 ;
  protonLowEnergy      = 1.*keV ; 
  protonHighEnergy     = 2.*MeV ;
  antiProtonLowEnergy  = 1.*keV ;
  antiProtonHighEnergy = 2.*MeV ;
  minGammaEnergy       = 25.*keV;
  minElectronEnergy    = 25.*keV;
  verboseLevel         = 0;
  shellCS = new G4hShellCrossSection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hLowEnergyIonisation::~G4hLowEnergyIonisation() 
{
  if (theMeanFreePathTable) {
    theMeanFreePathTable->clearAndDestroy();
    delete theMeanFreePathTable;
  }
  if(theBetheBlochModel)delete theBetheBlochModel;
  if(theProtonModel)delete theProtonModel;
  if(theAntiProtonModel)delete theAntiProtonModel;
  if(theNuclearStoppingModel)delete theNuclearStoppingModel;
  if(theIonEffChargeModel)delete theIonEffChargeModel;
  if(theIonChuFluctuationModel)delete theIonChuFluctuationModel;
  if(theIonYangFluctuationModel)delete theIonYangFluctuationModel;
  if(shellVacancy) delete shellVacancy;
  if(shellCS) delete shellCS;
  cutForDelta.clear();
  G4int length = zFluoDataVector.size();
  if(length) {
    for(G4int i=0; i<length; i++) {
      delete &(zFluoDataVector[i]);
    }
    zFluoDataVector.clear();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::SetElectronicStoppingPowerModel(
                             const G4ParticleDefinition* aParticle,
                             const G4String& dedxTable)
  // This method defines the ionisation parametrisation method via its name 
{
  if(0 < aParticle->GetPDGCharge()) {
    SetProtonElectronicStoppingPowerModel(dedxTable) ;
  } else {
    SetAntiProtonElectronicStoppingPowerModel(dedxTable) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::InitializeParametrisation() 

{
  G4Proton* theProton = G4Proton::Proton();
  G4AntiProton* theAntiProton = G4AntiProton::AntiProton();

  // Define models for parametrisation of electronic energy losses
  theBetheBlochModel = new G4hBetheBlochModel("Bethe-Bloch") ;
  theProtonModel = new G4hParametrisedLossModel(theProtonTable) ;
  theAntiProtonModel = new G4QAOLowEnergyLoss(theAntiProtonTable) ;
  theNuclearStoppingModel = new G4hNuclearStoppingModel(theNuclearTable) ;
  theIonEffChargeModel = new G4hIonEffChargeSquare("Ziegler1988") ;
  theIonChuFluctuationModel = new G4IonChuFluctuationModel("Chu") ;
  theIonYangFluctuationModel = new G4IonYangFluctuationModel("Yang") ;

  // Energy limits for parametrisation of electronic energy losses
  protonLowEnergy = G4std::max(protonLowEnergy,
                               theProtonModel->LowEnergyLimit(theProton)) ;   

  G4double x1 = theBetheBlochModel->LowEnergyLimit(theProton) ;   
  G4double x2 = theProtonModel->HighEnergyLimit(theProton) ;   
  if(protonHighEnergy < x1) protonHighEnergy = x1 ;
  if(protonHighEnergy > x2) protonHighEnergy = x2 ;
  
  antiProtonLowEnergy = G4std::max(antiProtonLowEnergy,
                      theAntiProtonModel->LowEnergyLimit(theAntiProton)) ;   

  x2 = theAntiProtonModel->HighEnergyLimit(theAntiProton) ;   
  if(antiProtonHighEnergy < x1) antiProtonHighEnergy = x1 ;
  if(antiProtonHighEnergy > x2) antiProtonHighEnergy = x2 ;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::BuildPhysicsTable(
                       const G4ParticleDefinition& aParticleType) 
  
  //  just call BuildLossTable+BuildLambdaTable
{
  if(verboseLevel > 0) {
    G4cout << "G4hLowEnergyIonisation::BuildPhysicsTable for "
           << aParticleType.GetParticleName()
           << " mass(MeV)= " << aParticleType.GetPDGMass()/MeV
           << " charge= " << aParticleType.GetPDGCharge()/eplus
           << " type= " << aParticleType.GetParticleType()
           << G4endl;

    if(verboseLevel > 1) {
      G4ProcessVector* pv = aParticleType.GetProcessManager()->GetProcessList();
     
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

  if(aParticleType.GetParticleType() == "nucleus" && 
     aParticleType.GetParticleName() != "GenericIon" &&
     theMeanFreePathTable) {

     G4EnergyLossTables::Register(&aParticleType,  
              theDEDXpTable,
              theRangepTable,
              theInverseRangepTable,
              theLabTimepTable,
              theProperTimepTable,
              LowestKineticEnergy, HighestKineticEnergy,
              proton_mass_c2/aParticleType.GetPDGMass(),
              TotBin);

     return;
  }
  
  InitializeParametrisation() ;
  G4Proton* theProton = G4Proton::Proton();
  G4AntiProton* theAntiProton = G4AntiProton::AntiProton();

  charge = aParticleType.GetPDGCharge()/eplus;
  chargeSquare = charge*charge ;

  // ---- MGP ---- workaround for the deprecated "cuts per material"
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  const G4Material* material = (*theMaterialTable)[0];  
  G4double electronCutInRange = G4Electron::Electron()->GetRangeThreshold(material);
  //                    was   = G4Electron::Electron()->GetCuts(); 
  // ---- MGP ----

  // Define cuts

  //  create table
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  cutForDelta.clear();
  cutForGamma.clear();

  for (G4int j=0; j<numOfMaterials; j++) {
      
    // get material parameters needed for the energy loss calculation  
    const G4Material* material= (*theMaterialTable)[j];

    // the cut cannot be below lowest limit
    G4double tCut = G4Electron::Electron()->GetEnergyThreshold(material);
    if(tCut > HighestKineticEnergy) tCut = HighestKineticEnergy;

    G4double excEnergy = material->GetIonisation()->GetMeanExcitationEnergy();

    tCut = G4std::max(tCut,excEnergy);
    cutForDelta.push_back(tCut);

    // the cut cannot be below lowest limit
    tCut = G4Gamma::Gamma()->GetEnergyThreshold(material);
    if(tCut > HighestKineticEnergy) tCut = HighestKineticEnergy;
    tCut = G4std::max(tCut,minGammaEnergy);
    cutForGamma.push_back(tCut);
  }

  if(verboseLevel > 0) {
    G4cout << "Cuts are defined " << G4endl;
  }


  if(0.0 < charge) 
    {

      if( (ptableElectronCutInRange != electronCutInRange)  
	  || (theDEDXpTable == 0))
	{
	  BuildLossTable(*theProton) ;
	  RecorderOfpProcess[CounterOfpProcess] = theLossTable ;
	  CounterOfpProcess++;
	}

    } else{

      if( (pbartableElectronCutInRange != electronCutInRange)  
	  || (theDEDXpbarTable == 0))
	{
	  BuildLossTable(*theAntiProton) ;
	  RecorderOfpbarProcess[CounterOfpbarProcess] = theLossTable ;
	  CounterOfpbarProcess++;
	}
    }

  if(verboseLevel > 0) {
    G4cout << "G4hLowEnergyIonisation::BuildPhysicsTable: "
           << "Loss table is built" << G4endl;
  }
  
  BuildLambdaTable(aParticleType) ;
  BuildDataForFluorescence(aParticleType);

  if(verboseLevel > 0) {
    G4cout << "G4hLowEnergyIonisation::BuildPhysicsTable: "
           << "DEDX table will be built" << G4endl;
  }

  BuildDEDXTable(aParticleType) ;

  if((&aParticleType == theProton) ) PrintInfoDefinition() ;

  if(verboseLevel > 0) {
    G4cout << "G4hLowEnergyIonisation::BuildPhysicsTable: end for "
           << aParticleType.GetParticleName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::BuildLossTable(
                             const G4ParticleDefinition& aParticleType) 
{
  
  // Initialisation
  G4double lowEdgeEnergy , ionloss, ionlossBB, paramB ;
  G4double lowEnergy, highEnergy;
  G4Proton* theProton = G4Proton::Proton();

  if(aParticleType == *theProton) {
    lowEnergy = protonLowEnergy ;
    highEnergy = protonHighEnergy ;
    charge = 1.0 ;
  } else {
    lowEnergy = antiProtonLowEnergy ;
    highEnergy = antiProtonHighEnergy ;
    charge = -1.0 ;
  }
  chargeSquare = 1.0 ;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  
  //  create table
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();
  
  if ( theLossTable) {
    theLossTable->clearAndDestroy();
    delete theLossTable;
  }

  theLossTable = new G4PhysicsTable(numOfMaterials);
  
  //  loop for materials  
  for (G4int j=0; j<numOfMaterials; j++) {

    // create physics vector and fill it  
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(LowestKineticEnergy, 
                                                         HighestKineticEnergy, 
                                                         TotBin);
      
    // get material parameters needed for the energy loss calculation  
    const G4Material* material= (*theMaterialTable)[j];

    // low energy of Bethe-Bloch formula for this material
    G4double highE = G4std::max(highEnergy,theBetheBlochModel->
                                LowEnergyLimit(&aParticleType,material)) ;

    if ( charge > 0.0 ) {
      ionloss = ProtonParametrisedDEDX(material,highE) ;
    } else {
      ionloss = AntiProtonParametrisedDEDX(material,highE) ;
    }

    ionlossBB = theBetheBlochModel->TheValue(&aParticleType,material,highE) ;
    ionlossBB -= DeltaRaysEnergy(material,highE,proton_mass_c2) ;

    /*
    if(theBarkas) {
      ionlossBB += BarkasTerm(material,highE)*charge ;
      ionlossBB += BlochTerm(material,highE,1.0) ;
    }
    */
    
    paramB =  ionloss/ionlossBB - 1.0 ; 

    // now comes the loop for the kinetic energy values      
    for (G4int i = 0 ; i < TotBin ; i++) {
      lowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;

      // low energy part for this material, parametrised energy loss formulae
      if ( lowEdgeEnergy < highE ) {
	 
        if ( charge > 0.0 ) {
          ionloss = ProtonParametrisedDEDX(material,lowEdgeEnergy) ;
	} else {
          ionloss = AntiProtonParametrisedDEDX(material,lowEdgeEnergy) ;
	}

      } else {
	    
        // high energy part for this material, Bethe-Bloch formula
        ionloss = theBetheBlochModel->TheValue(theProton,material, 
                                               lowEdgeEnergy) ; 

        ionloss -= DeltaRaysEnergy(material,lowEdgeEnergy,proton_mass_c2) ;

	/*
        if(theBarkas) {
          ionloss += BarkasTerm(material,lowEdgeEnergy)*charge ;
          ionloss += BlochTerm(material,lowEdgeEnergy,1.0) ;
        }
	*/
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::BuildDataForFluorescence(
                         const G4ParticleDefinition& aParticleType)
{

  if(verboseLevel > 1) {
    G4cout << "G4hLowEnergyIonisation::BuildDataForFluorescence for "
           << aParticleType.GetParticleName() << " is started" << G4endl;
  }

  // fill data for fluorescence

  G4double mass = aParticleType.GetPDGMass();
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  
  //  create table
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  if (shellVacancy != 0) delete shellVacancy;
  shellVacancy = new G4ShellVacancy();
  G4DataVector* ksi = 0;
  G4DataVector* ksi1 = 0;
  G4DataVector* energy = 0;
  G4DataVector* energy1 = 0;
  size_t binForFluo = TotBin/10;
  G4int length = zFluoDataVector.size();
  if(length > 0) {
    for(G4int i=0; i<length; i++) {
      G4VEMDataSet* x = zFluoDataVector[i];
      delete x;
    }
    zFluoDataVector.clear();
  }

  G4PhysicsLogVector* bVector = new G4PhysicsLogVector(LowestKineticEnergy,
		                		       HighestKineticEnergy,
						       binForFluo);
  const G4AtomicTransitionManager* transitionManager = 
                             G4AtomicTransitionManager::Instance();

  G4double bindingEnergy;
  //  G4double x;
  //  G4double y;

  //  loop for materials  
  for (G4int j=0; j<numOfMaterials; j++) {

    // get material parameters needed for the energy loss calculation  
    const G4Material* material= (*theMaterialTable)[j];

    const G4ElementVector* theElementVector = material->GetElementVector();
    size_t NumberOfElements = material->GetNumberOfElements() ;
    const G4double* theAtomicNumDensityVector = 
                    material->GetAtomicNumDensityVector();
    G4VDataSetAlgorithm* interp = new G4SemiLogInterpolation();
    G4VEMDataSet* xsis = new G4CompositeEMDataSet(interp, 1., 1.);
    G4VDataSetAlgorithm* interp1 = new G4SemiLogInterpolation();
    G4VEMDataSet* xsis1 = new G4CompositeEMDataSet(interp1, 1., 1.);
    
    G4double tCut = cutForDelta[j];
    G4double elDensity = 1.;

    for (size_t iel=0; iel<NumberOfElements; iel++ ) {

      G4int Z = (G4int)((*theElementVector)[iel]->GetZ());
      G4int nShells = transitionManager->NumberOfShells(Z);
      energy = new G4DataVector();
      ksi    = new G4DataVector();
      energy1= new G4DataVector();
      ksi1   = new G4DataVector();
      //if(NumberOfElements > 1) 
      elDensity = theAtomicNumDensityVector[iel]/((G4double)nShells);

      for (size_t j = 0; j<binForFluo; j++) {

        G4double tkin  = bVector->GetLowEdgeEnergy(j);
        G4double gamma = tkin/mass + 1.;
        G4double beta2 = 1.0 - 1.0/(gamma*gamma);
        G4double r     = electron_mass_c2/mass;
        G4double tmax  = 2.*electron_mass_c2*(gamma*gamma - 1.)/(1. + 2.*gamma*r + r*r);
        G4double cross   = 0.;
        G4double cross1  = 0.;
        G4double eAverage= 0.;
        G4double tmin = G4std::min(tCut,tmax);
        G4double rel;

        for (G4int n=0; n<nShells; n++) {

          bindingEnergy = transitionManager->Shell(Z, n)->BindingEnergy();
          if (tmin > bindingEnergy) {
            rel = log(tmin/bindingEnergy);
            eAverage   += rel - beta2*(tmin - bindingEnergy)/tmax;
            cross      += 1.0/bindingEnergy - 1.0/tmin - beta2*rel/tmax;
	  }
          if (tmax > tmin) {
	    cross1     += 1.0/tmin - 1.0/tmax - beta2*log(tmax/tmin)/tmax;
	  }
	}

        cross1 *= elDensity;
        energy1->push_back(tkin);
        ksi1->push_back(cross1);

        if(eAverage > 0.) cross /= eAverage;
        else              cross  = 0.;

        energy->push_back(tkin);
        ksi->push_back(cross);
      }
      G4VDataSetAlgorithm* algo = interp->Clone();
      G4VEMDataSet* set = new G4EMDataSet(Z,energy,ksi,algo,1.,1.);
      xsis->AddComponent(set);
      G4VDataSetAlgorithm* algo1 = interp1->Clone();
      G4VEMDataSet* set1 = new G4EMDataSet(Z,energy1,ksi1,algo1,1.,1.);
      xsis1->AddComponent(set1);
    }
    if(verboseLevel > 1) {
      G4cout << "### Shell inverse cross sections for " 
             << material->GetName() << G4endl;
      xsis->PrintData();
      G4cout << "### Atom cross sections for " 
             << material->GetName() << G4endl;
      xsis1->PrintData();
    }
    shellVacancy->AddXsiTable(xsis);
    zFluoDataVector.push_back(xsis1);    
  }
  delete bVector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::BuildLambdaTable(
                       const G4ParticleDefinition& aParticleType) 
  
{
  // Build mean free path tables for the delta ray production process
  //     tables are built for MATERIALS 
  
  if(verboseLevel > 1) {
    G4cout << "G4hLowEnergyIonisation::BuildLambdaTable for "
           << aParticleType.GetParticleName() << " is started" << G4endl;
  }


  G4double lowEdgeEnergy, value;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  charge = aParticleType.GetPDGCharge()/eplus ;
  chargeSquare = charge*charge ;
  initialMass = aParticleType.GetPDGMass();
    
  //create table
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();
  
  if (theMeanFreePathTable) {
    theMeanFreePathTable->clearAndDestroy();
    delete theMeanFreePathTable;
  }
  
  theMeanFreePathTable = new G4PhysicsTable(numOfMaterials);
    
  // loop for materials 
  
  for (G4int J=0 ; J < numOfMaterials; J++) { 

    //create physics vector then fill it ....      
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(LowestKineticEnergy, 
                                                         HighestKineticEnergy, 
                                                         TotBin);
      
    // compute the (macroscopic) cross section first
      
    const G4Material* material = (*theMaterialTable)[J];   
    const G4ElementVector* theElementVector = 
                           material->GetElementVector() ;
    const G4double* theAtomicNumDensityVector = 
                           material->GetAtomicNumDensityVector();
    const G4int NumberOfElements = material->GetNumberOfElements() ;
      
      // get the electron kinetic energy cut for the actual material,
      //  it will be used in ComputeMicroscopicCrossSection
      // ( it is the SAME for ALL the ELEMENTS in THIS MATERIAL )
      //   ------------------------------------------------------
      
    G4double deltaCut = cutForDelta[J];
      
    for ( G4int i = 0 ; i < TotBin ; i++ ) {
      lowEdgeEnergy = aVector->GetLowEdgeEnergy(i) ;
      G4double sigma = 0.0 ;  
      for (G4int iel=0; iel<NumberOfElements; iel++ ) {
	sigma += theAtomicNumDensityVector[iel]*
		 ComputeMicroscopicCrossSection(
                        aParticleType,
			lowEdgeEnergy,
                        (*theElementVector)[iel]->GetZ(),
                        deltaCut ) ;
      }
	  
      // mean free path = 1./macroscopic cross section
	  
      value = sigma<=0 ? DBL_MAX : 1./sigma ;     
	  
      aVector->PutValue(i, value) ;
    }
  
    theMeanFreePathTable->insert(aVector);
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::ComputeMicroscopicCrossSection(
                           const G4ParticleDefinition& aParticleType,
			         G4double kineticEnergy,
				 G4double atomicNumber,
				 G4double deltaCutInEnergy) const
{
  //******************************************************************
  // cross section formula is OK for spin=0, 1/2, 1 only !
  // *****************************************************************
  
  // calculates the microscopic cross section in GEANT4 internal units
  //    ( it is called for elements , AtomicNumber = z )
  
  G4double energy, gamma, beta2, tmax, var;
  G4double totalCrossSection = 0.0 ;
  
  G4double particleMass = initialMass;
  
  // get particle data ...................................
  
  energy = kineticEnergy + particleMass;
  
  // some kinematics......................
  
  gamma  = energy/particleMass;
  beta2  = 1.0 - 1.0/(gamma*gamma);
  var    = electron_mass_c2/particleMass;
  tmax   = 2.*electron_mass_c2*(gamma*gamma - 1.)/(1. + 2.*gamma*var + var*var);
  
  // now you can calculate the total cross section 
  
  if( tmax > deltaCutInEnergy ) {
    
    var=deltaCutInEnergy/tmax;
    totalCrossSection = (1.0 - var*(1.0 - beta2*log(var))) / deltaCutInEnergy ;
    G4double spin = aParticleType.GetPDGSpin() ;
    
    // +term for spin=1/2 particle
    if( 0.5 == spin )
      totalCrossSection +=  0.5 * (tmax - deltaCutInEnergy) / (energy*energy);
    
    // +term for spin=1 particle
    else if( 0.9 < spin )
      totalCrossSection += -log(var)/(3.0*deltaCutInEnergy) +
	(tmax - deltaCutInEnergy) * ( (5.0+ 1.0/var)*0.25 / (energy*energy) - 
	 beta2 / (tmax * deltaCutInEnergy) ) / 3.0 ;
    
    totalCrossSection *= twopi_mc2_rcl2 * atomicNumber / beta2 ;
  }
  
  return totalCrossSection ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetMeanFreePath(const G4Track& trackData,
                                                   G4double previousStepSize,
                                              enum G4ForceCondition* condition)
{
   const G4DynamicParticle* aParticle = trackData.GetDynamicParticle() ;
   G4Material* aMaterial = trackData.GetMaterial() ;
   G4double meanFreePath;
   G4bool isOutRange ;
 
   *condition = NotForced ;

   G4double kineticEnergy = (aParticle->GetKineticEnergy())*initialMass/(aParticle->GetMass());
   charge = aParticle->GetCharge();
   chargeSquare = theIonEffChargeModel->TheValue(aParticle, aMaterial);

   if(kineticEnergy < LowestKineticEnergy) meanFreePath = DBL_MAX;

   else {
     if(kineticEnergy > HighestKineticEnergy)
                    kineticEnergy = HighestKineticEnergy;
     meanFreePath = (((*theMeanFreePathTable)(aMaterial->GetIndex()))->
                    GetValue(kineticEnergy,isOutRange))/chargeSquare;
     }

   return meanFreePath ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::GetConstraints(
                                 const G4DynamicParticle* particle,
				 const G4Material* material) 
{
  // returns the Step limit
  // dEdx is calculated as well as the range  
  // based on Effective Charge Approach
 
  G4Proton* theProton = G4Proton::Proton();
  G4AntiProton* theAntiProton = G4AntiProton::AntiProton();

  G4double stepLimit = 0.0 ;
  G4double dx, highEnergy;
  
  G4double massRatio = proton_mass_c2/(particle->GetMass()) ;
  G4double kineticEnergy = particle->GetKineticEnergy() ;
  
  // Scale the kinetic energy
  
  G4double tscaled = kineticEnergy*massRatio ; 
  
  if(charge > 0.0) {
    
    highEnergy = protonHighEnergy ;

    fRangeNow = G4EnergyLossTables::GetRange(theProton, tscaled, material);
    dx = G4EnergyLossTables::GetRange(theProton, highEnergy, material);
    fdEdx = G4EnergyLossTables::GetDEDX(theProton, tscaled, material) 
          * chargeSquare ;

    if(tscaled > highEnergy) {
        // Correction for positive ions
      if(theBarkas) {
        fdEdx += BarkasTerm(material,tscaled)*sqrt(chargeSquare)*chargeSquare; 
        fdEdx += BlochTerm(material,tscaled,chargeSquare); 
      }
    }
    
    // Antiprotons and negative hadrons 
  } else {

    highEnergy = antiProtonHighEnergy ;
    fRangeNow = G4EnergyLossTables::GetRange(theAntiProton, tscaled, material);
    dx = G4EnergyLossTables::GetRange(theAntiProton, highEnergy, material);
    fdEdx = G4EnergyLossTables::GetDEDX(theAntiProton, tscaled, material) 
          * chargeSquare ;

    if(tscaled > highEnergy) {

      // Correction for positive ions
      if(theBarkas) {
        fdEdx -= BarkasTerm(material,tscaled)*sqrt(chargeSquare)*chargeSquare; 
        fdEdx += BlochTerm(material,tscaled,chargeSquare); 
      }
    }
  }
  
  // scaling back
  fRangeNow /= (chargeSquare*massRatio) ;
  dx        /= (chargeSquare*massRatio) ;
  stepLimit  = fRangeNow ;
  
  // compute the (random) Step limit in standard energy range
  if(tscaled > highEnergy ) {
    if(fRangeNow > finalRange) {
      stepLimit = (c1lim*fRangeNow+c2lim+c3lim/fRangeNow) ;
      
      //  randomise this value
      if(rndmStepFlag) stepLimit = finalRange + 
                                  (stepLimit-finalRange)*G4UniformRand() ;
    }
    // cut step in front of Bragg's peak
    if(stepLimit > fRangeNow - dx*0.9) stepLimit = fRangeNow - dx*0.9 ;

  // Step limit in low energy range
  } else {
    stepLimit = G4std::min(fRangeNow, dx*paramStepLimit) ;    
  }

  return stepLimit ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4hLowEnergyIonisation::AlongStepDoIt(
                                           const G4Track& trackData, 
                                           const G4Step& stepData) 
{
  // compute the energy loss after a step 
  G4Proton* theProton = G4Proton::Proton();
  G4AntiProton* theAntiProton = G4AntiProton::AntiProton();
  G4double finalT = 0.0 ;
  
  aParticleChange.Initialize(trackData) ;

  G4Material* material = trackData.GetMaterial() ;
  
  // get the actual (true) Step length from stepData 
  const G4double step = stepData.GetStepLength() ;
  
  const G4DynamicParticle* particle = trackData.GetDynamicParticle() ;
  
  G4double kineticEnergy = particle->GetKineticEnergy() ;
  G4double massRatio = proton_mass_c2/(particle->GetMass()) ;
  G4double tscaled= kineticEnergy*massRatio ; 
  G4double eloss = 0.0 ;
  G4double nloss = 0.0 ;

  
    // very small particle energy
  if(kineticEnergy < MinKineticEnergy) {

    eloss = kineticEnergy ;  

    // particle energy outside tabulated energy range
  } else if( kineticEnergy > HighestKineticEnergy) { 
    eloss = step*fdEdx ; 
 
    // big step
  } else if(step >= fRangeNow ) {
    eloss = kineticEnergy ;

    // tabulated range 
  } else {
    
    // step longer than linear step limit 
    if(step > linLossLimit*fRangeNow) {
      
      G4double rscaled= fRangeNow*massRatio*chargeSquare ;
      G4double sscaled=   step   *massRatio*chargeSquare ;
      
      if(charge > 0.0) {
        eloss = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                    theProton,rscaled, material) -
	        G4EnergyLossTables::GetPreciseEnergyFromRange(
                                    theProton,rscaled-sscaled,material) ;
       
      } else {
        eloss = G4EnergyLossTables::GetPreciseEnergyFromRange(
                                    theAntiProton,rscaled,material) -
	        G4EnergyLossTables::GetPreciseEnergyFromRange(
                                    theAntiProton,rscaled-sscaled,material) ;
      }
      eloss /= massRatio ;

    // step shorter than linear step limit       
    } else {
      eloss = step*fdEdx ;
    }
    // Correction for positive ions
    if(theBarkas && 1.0 < charge) {
      G4double ts = tscaled - eloss*0.5*massRatio;
      if(ts < protonHighEnergy) ts = protonHighEnergy;
      eloss += BarkasTerm(material,ts)*charge*chargeSquare*step;
      eloss += BlochTerm(material,ts,chargeSquare)*step; 
    }
    if(nStopping && tscaled < protonHighEnergy) {
      nloss = (theNuclearStoppingModel->TheValue(particle, material))*step;
    }
  }

  if(eloss < 0.0) eloss = 0.0;
  
  finalT = kineticEnergy - eloss - nloss;

  if( EnlossFlucFlag && 0.0 < eloss && finalT > MinKineticEnergy) {
    
    //  now the electron loss with fluctuation
    eloss = ElectronicLossFluctuation(particle, material, eloss, step) ;
    if(eloss < 0.0) eloss = 0.0;
    finalT = kineticEnergy - eloss - nloss;    
  }
  
  //  stop particle if the kinetic energy <= MinKineticEnergy
  if (finalT <= MinKineticEnergy ) {
     
     finalT = 0.0;
      if(!particle->GetDefinition()->GetProcessManager()->
                     GetAtRestProcessVector()->size())
        aParticleChange.SetStatusChange(fStopAndKill);
      else 
        aParticleChange.SetStatusChange(fStopButAlive);
  } 

  aParticleChange.SetEnergyChange( finalT );
  G4double edep = kineticEnergy-finalT;

  // Deexcitation only of ionised atoms
  eloss = G4std::min(edep, eloss);
  
  G4double hMass = particle->GetMass();
  G4std::vector<G4DynamicParticle*>* newpart = 0;
  G4DynamicParticle* part = 0;

  if(theFluo) newpart = DeexciteAtom(material, kineticEnergy, hMass, eloss);

  if(newpart != 0) {

    size_t nSecondaries = newpart->size();
    aParticleChange.SetNumberOfSecondaries(nSecondaries);
    G4Track* newtrack = 0;
    const G4StepPoint* preStep = stepData.GetPreStepPoint();
    const G4StepPoint* postStep = stepData.GetPostStepPoint();
    G4ThreeVector r = preStep->GetPosition();
    G4ThreeVector deltaR = postStep->GetPosition();
    deltaR -= r;
    G4double t = preStep->GetGlobalTime();
    G4double deltaT = postStep->GetGlobalTime();
    deltaT -= t;
    G4double time, q, e;
    G4ThreeVector position;

    for(size_t i=0; i<nSecondaries; i++) {

      part = (*newpart)[i]; 
      if(part) {

        e = part->GetKineticEnergy();
        if(e <= edep) {

          edep -= e;
          q = G4UniformRand();
          time = deltaT*q + t;
          position  = deltaR*q;
          position += r;
          newtrack = new G4Track(part, time, position);
          aParticleChange.AddSecondary(newtrack);
	
        } else {
         
          delete part;

	}
      }
    }
    delete newpart;   
  }

  aParticleChange.SetLocalEnergyDeposit(edep);
  return &aParticleChange ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::ProtonParametrisedDEDX(
                                 const G4Material* material, 
				       G4double kineticEnergy) const
{
  G4Proton* theProton = G4Proton::Proton();
  G4double eloss = 0.0;

    // Free Electron Gas Model  
  if(kineticEnergy < protonLowEnergy) {
    eloss = (theProtonModel->TheValue(theProton, material, protonLowEnergy))
          * sqrt(kineticEnergy/protonLowEnergy) ;
    
    // Parametrisation
  } else {
    eloss = theProtonModel->TheValue(theProton, material, kineticEnergy) ;
  }
  
  // Delta rays energy
  eloss -= DeltaRaysEnergy(material,kineticEnergy,proton_mass_c2) ;

  if(verboseLevel > 2) {
    G4cout << "p E(MeV)= " << kineticEnergy/MeV
           << " dE/dx(MeV/mm)= " << eloss*mm/MeV 
           << " for " << material->GetName() 
           << " model: " << theProtonModel << G4endl;
  }

  if(eloss < 0.0) eloss = 0.0 ;  
  
  return eloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::AntiProtonParametrisedDEDX(
                                 const G4Material* material, 
				       G4double kineticEnergy) const
{
  G4AntiProton* theAntiProton = G4AntiProton::AntiProton();
  G4double eloss = 0.0 ;

  // Antiproton model is used
  if(theAntiProtonModel->IsInCharge(theAntiProton,material)) {
    if(kineticEnergy < antiProtonLowEnergy) {
      eloss = theAntiProtonModel->TheValue(theAntiProton,material,
                                           antiProtonLowEnergy);
    
    // Parametrisation 
    } else {
      eloss = theAntiProtonModel->TheValue(theAntiProton,material,
                                           kineticEnergy);
    }

  // The proton model is used + Barkas correction
  } else { 
    if(kineticEnergy < antiProtonLowEnergy) {
      eloss = theProtonModel->TheValue(G4Proton::Proton(),material,
                                       antiProtonLowEnergy);
    
    // Parametrisation 
    } else {
      eloss = theProtonModel->TheValue(G4Proton::Proton(),material,
                                       kineticEnergy);
    }
    if(theBarkas) eloss -= 2.0*BarkasTerm(material, kineticEnergy);
  }
  
  // Delta rays energy
  eloss -= DeltaRaysEnergy(material,kineticEnergy,proton_mass_c2) ;

  if(verboseLevel > 2) {
    G4cout << "pbar E(MeV)= " << kineticEnergy/MeV
           << " dE/dx(MeV/mm)= " << eloss*mm/MeV 
           << " for " << material->GetName() 
           << " model: " << theProtonModel << G4endl;
  }

  if(eloss < 0.0) eloss = 0.0 ;  
  
  return eloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::DeltaRaysEnergy(
                                 const G4Material* material, 
                                       G4double kineticEnergy,
                                       G4double particleMass) const
{
  G4double dloss = 0.0 ;

  G4double deltaCutNow = cutForDelta[(material->GetIndex())] ;   
  G4double electronDensity = material->GetElectronDensity();
  G4double eexc = material->GetIonisation()->GetMeanExcitationEnergy();

  G4double tau = kineticEnergy/particleMass ;    
  G4double rateMass = electron_mass_c2/particleMass ;
  
  // some local variables 
  
  G4double gamma,bg2,beta2,tmax,x ;
  
  gamma = tau + 1.0 ;
  bg2 = tau*(tau+2.0) ;
  beta2 = bg2/(gamma*gamma) ;
  tmax = 2.*electron_mass_c2*bg2/(1.0+2.0*gamma*rateMass+rateMass*rateMass) ;

  // Validity range for delta electron cross section
  G4double deltaCut = G4std::max(deltaCutNow, eexc);
    
  if ( deltaCut < tmax) {
    x = deltaCut / tmax ;
    dloss = ( beta2 * (x - 1.0) - log(x) ) * twopi_mc2_rcl2 
          * electronDensity / beta2 ;
  }
  return dloss ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4hLowEnergyIonisation::PostStepDoIt(
                                           const G4Track& trackData,   
					   const G4Step& stepData) 
{
  // Units are expressed in GEANT4 internal units.
  
  G4double KineticEnergy,TotalEnergy,TotalMomentum,betasquare,
           DeltaKineticEnergy,DeltaTotalMomentum,costheta,sintheta,phi,
           dirx,diry,dirz,finalKineticEnergy,finalPx,finalPy,finalPz,
           x,xc,grej,Psquare,Esquare,rate,finalMomentum ;
  
  aParticleChange.Initialize(trackData) ;
  G4Material* aMaterial = trackData.GetMaterial() ;
  
  const G4DynamicParticle* aParticle = trackData.GetDynamicParticle() ;
  
  // some kinematics

  ParticleMass=aParticle->GetDefinition()->GetPDGMass();
  KineticEnergy=aParticle->GetKineticEnergy();
  TotalEnergy=KineticEnergy + ParticleMass ;
  Psquare=KineticEnergy*(TotalEnergy+ParticleMass) ;
  Esquare=TotalEnergy*TotalEnergy;
  betasquare=Psquare/Esquare;
  G4ThreeVector ParticleDirection = aParticle->GetMomentumDirection() ;
  
  G4double gamma= KineticEnergy/ParticleMass + 1.;
  G4double r    = electron_mass_c2/ParticleMass;
  G4double tmax = 2.*electron_mass_c2*(gamma*gamma - 1.)/(1. + 2.*gamma*r + r*r);
  
  // Validity range for delta electron cross section
  G4double DeltaCut = cutForDelta[aMaterial->GetIndex()];

  // This should not be a case
  if(DeltaCut >= tmax)
       return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
    
  xc   = DeltaCut / tmax;
  rate = tmax / TotalEnergy;
  rate = rate*rate ;
  G4double spin = aParticle->GetDefinition()->GetPDGSpin() ;   
      
  // sampling follows ...      
      do {
	x=xc/(1.-(1.-xc)*G4UniformRand());

        if(0.0 == spin) {
          grej = 1.0 - betasquare * x ;

        } else if (0.5 == spin) {
          grej = (1.0 - betasquare * x + 0.5*x*x*rate) / (1.0 + 0.5 * rate) ;

        } else {
          grej = (1.0 - betasquare * x ) * (1.0 + x/ (3.0*xc)) + 
	    x * x * rate * (1.0 + 0.5 * x / xc) / 3.0 / 
	    (1.0 + 1.0/(3.0*xc) + rate *(1.0+ 0.5/xc) /3.0) ;
        }

      } while( G4UniformRand() > grej );
  
  
  DeltaKineticEnergy = x * tmax;
    
  DeltaTotalMomentum = sqrt(DeltaKineticEnergy * (DeltaKineticEnergy +
						  2. * electron_mass_c2 )) ;
  TotalMomentum = sqrt(Psquare) ;
  costheta = DeltaKineticEnergy * (TotalEnergy + electron_mass_c2)
    /(DeltaTotalMomentum * TotalMomentum) ;
  
  //  protection against costheta > 1 or < -1   ---------------
  if ( costheta < -1. ) 
    costheta = -1. ;
  if ( costheta > +1. ) 
    costheta = +1. ;
  
  //  direction of the delta electron  ........
  phi = twopi * G4UniformRand() ; 
  sintheta = sqrt(1. - costheta*costheta);
  dirx = sintheta * cos(phi) ;
  diry = sintheta * sin(phi) ;
  dirz = costheta ;
  
  G4ThreeVector DeltaDirection(dirx,diry,dirz) ;
  DeltaDirection.rotateUz(ParticleDirection) ;
  
  // create G4DynamicParticle object for delta ray
  G4DynamicParticle *theDeltaRay = new G4DynamicParticle;
  theDeltaRay->SetKineticEnergy( DeltaKineticEnergy );
  theDeltaRay->SetMomentumDirection(DeltaDirection.x(),
                                    DeltaDirection.y(),
                                    DeltaDirection.z()); 
  theDeltaRay->SetDefinition(G4Electron::Electron());
  
  // fill aParticleChange 
  finalKineticEnergy = KineticEnergy - DeltaKineticEnergy ;

  // Generation of Fluorescence and Auger
  size_t nSecondaries = 0;
  size_t totalNumber  = 1;
  G4std::vector<G4DynamicParticle*>* secondaryVector = 0;
  G4DynamicParticle* aSecondary = 0;
  G4ParticleDefinition* type = 0;

  // Select atom and shell
  G4int Z = SelectRandomAtom(aMaterial, KineticEnergy);

  if(theFluo && Z > 5) {
    G4int shell = shellCS->SelectRandomShell(Z, KineticEnergy,
                                             ParticleMass,DeltaKineticEnergy);
    const G4AtomicShell* atomicShell = 
                (G4AtomicTransitionManager::Instance())->Shell(Z, shell);
    G4double bindingEnergy = atomicShell->BindingEnergy();

    if(verboseLevel > 1) {
      G4cout << "PostStep Z= " << Z << " shell= " << shell
             << " bindingE(keV)= " << bindingEnergy/keV
             << " finalE(keV)= " << finalKineticEnergy/keV
             << G4endl;
    } 

    // Fluorescence data start from element 6
  
    if (finalKineticEnergy >= bindingEnergy 
         && (bindingEnergy >= minGammaEnergy 
         ||  bindingEnergy >= minElectronEnergy) ) {

      G4int shellId = atomicShell->ShellId();
      secondaryVector = deexcitationManager.GenerateParticles(Z, shellId);

      if (secondaryVector != 0) {

        nSecondaries = secondaryVector->size();
        for (size_t i = 0; i<nSecondaries; i++) {

          aSecondary = (*secondaryVector)[i];
          if (aSecondary) {
	  
            G4double e = aSecondary->GetKineticEnergy();
            type = aSecondary->GetDefinition();
            if (e < finalKineticEnergy && 
                 ((type == G4Gamma::Gamma() && e > minGammaEnergy ) || 
                  (type == G4Electron::Electron() && e > minElectronEnergy ))) {

              finalKineticEnergy -= e;
              totalNumber++;

	    } else {

              delete aSecondary;
              (*secondaryVector)[i] = 0;
	    }
	  }
	}
      }
    }
  }    
  
  // Save delta-electrons
 
  G4double edep = 0.0;

  if (finalKineticEnergy > MinKineticEnergy)
    {
      finalPx = TotalMomentum*ParticleDirection.x()
	- DeltaTotalMomentum*DeltaDirection.x();
      finalPy = TotalMomentum*ParticleDirection.y()
	- DeltaTotalMomentum*DeltaDirection.y();
      finalPz = TotalMomentum*ParticleDirection.z()
	- DeltaTotalMomentum*DeltaDirection.z();
      finalMomentum =
	sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz) ;
      finalPx /= finalMomentum ;
      finalPy /= finalMomentum ;
      finalPz /= finalMomentum ;
      
      aParticleChange.SetMomentumChange( finalPx,finalPy,finalPz );
    }
  else
    {
      edep = finalKineticEnergy;
      finalKineticEnergy = 0.;
      aParticleChange.SetMomentumChange(ParticleDirection.x(),
                      ParticleDirection.y(),ParticleDirection.z());
      if(!aParticle->GetDefinition()->GetProcessManager()->
                     GetAtRestProcessVector()->size())
        aParticleChange.SetStatusChange(fStopAndKill);
      else 
        aParticleChange.SetStatusChange(fStopButAlive);
    }
  
  aParticleChange.SetEnergyChange( finalKineticEnergy );
  aParticleChange.SetLocalEnergyDeposit (edep);
  aParticleChange.SetNumberOfSecondaries(totalNumber);
  aParticleChange.AddSecondary(theDeltaRay);

  // Save Fluorescence and Auger
  
  if (secondaryVector) {

    for (size_t l = 0; l < nSecondaries; l++) {

      aSecondary = (*secondaryVector)[l];
      if(aSecondary) {
        aParticleChange.AddSecondary(aSecondary);
      } 
    }
    delete secondaryVector;
  }
   
  return G4VContinuousDiscreteProcess::PostStepDoIt(trackData,stepData);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4std::vector<G4DynamicParticle*>* 
G4hLowEnergyIonisation::DeexciteAtom(const G4Material* material,
				           G4double incidentEnergy,
				           G4double hMass,
				           G4double eLoss)
{

  if (verboseLevel > 1) {
	G4cout << "DeexciteAtom: cutForPhotons(keV)= " << minGammaEnergy/keV 
               << "  cutForElectrons(keV)= " << minElectronEnergy/keV 
               << "  eLoss(MeV)= " << eLoss
               << G4endl;
  }

  if(eLoss < minGammaEnergy && eLoss < minElectronEnergy) return 0;

  G4int index    = material->GetIndex();
  //  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double gamma = incidentEnergy/hMass + 1;
  G4double beta2 = 1.0 - 1.0/(gamma*gamma);
  G4double r     = electron_mass_c2/hMass;
  G4double tmax  = 2.*electron_mass_c2*(gamma*gamma - 1.)/(1. + 2.*gamma*r + r*r);
  G4double tcut  = G4std::min(tmax,cutForDelta[index]);
  const G4AtomicTransitionManager* transitionManager = 
                             G4AtomicTransitionManager::Instance();

  size_t nElements = material->GetNumberOfElements();
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4bool stop = true;

  for (size_t j=0; j<nElements; j++) {

    G4int Z = (G4int)((*theElementVector)[j]->GetZ());
    G4double maxE = transitionManager->Shell(Z, 0)->BindingEnergy();

    if (Z > 5 && maxE < tcut && (maxE > minGammaEnergy || maxE > minElectronEnergy) ) {
      stop = false;
      break;
    }
  }

  if(stop) return 0;

  // create vector of tracks of secondary particles

  G4std::vector<G4DynamicParticle*>* partVector = 
         new G4std::vector<G4DynamicParticle*>;
  G4std::vector<G4DynamicParticle*>* secVector = 0;
  G4DynamicParticle* aSecondary = 0;
  G4ParticleDefinition* type = 0;
  G4double e, tkin, grej;
  G4ThreeVector position;
  G4int shell, shellId;

  // sample secondaries

  G4double etot = 0.0; 
  G4std::vector<G4int> n = shellVacancy->GenerateNumberOfIonisations(material,
                                         incidentEnergy, eLoss);

  for (size_t i=0; i<nElements; i++) {

    size_t nVacancies = n[i];
    G4int Z = (G4int)((*theElementVector)[i]->GetZ());
    G4double maxE = transitionManager->Shell(Z, 0)->BindingEnergy();

    if (nVacancies && Z  > 5 && maxE < tcut && (maxE > minGammaEnergy || maxE > minElectronEnergy)) {
      for(size_t j=0; j<nVacancies; j++) {
     
        // sampling follows       
        do {
	  tkin = tcut/(1.0 + (tcut/maxE - 1.0)*G4UniformRand());
          grej = 1.0 - beta2 * tkin/tmax;
  
        } while( G4UniformRand() > grej );

        shell = shellCS->SelectRandomShell(Z,incidentEnergy,hMass,tkin);

        shellId = transitionManager->Shell(Z, shell)->ShellId();
        G4double maxE = transitionManager->Shell(Z, shell)->BindingEnergy();
 
        if (maxE>minGammaEnergy || maxE>minElectronEnergy ) {
          secVector = deexcitationManager.GenerateParticles(Z, shellId);
	} else {
          secVector = 0;
	}

        if (secVector) {	

          for (size_t l = 0; l<secVector->size(); l++) {

            aSecondary = (*secVector)[l];
            if(aSecondary) {
	  
              e = aSecondary->GetKineticEnergy();
              type = aSecondary->GetDefinition();
              if ( etot + e <= eLoss &&
                   (type == G4Gamma::Gamma() && e > minGammaEnergy ) || 
                   (type == G4Electron::Electron() && e > minElectronEnergy)) {
                   
                     etot += e;  
                     partVector->push_back(aSecondary);

	      } else {
                     delete aSecondary;
	      }
	    }
	  }
          delete secVector;
	}
      }
    }
  }

  if(partVector->empty()) {
    delete partVector;
    return 0;
  }

  return partVector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4hLowEnergyIonisation::SelectRandomAtom(const G4Material* material, 
                                                     G4double kineticEnergy) const
{
  G4int nElements = material->GetNumberOfElements();
  G4int Z = 0;

  if(nElements == 1) {
    Z = (G4int)(material->GetZ());
    return Z;
  }

  const G4ElementVector* theElementVector = material->GetElementVector();
  G4std::vector<G4double> p;
  G4int index = material->GetIndex();

  G4double norm = 0.0;
  for (G4int j=0; j<nElements; j++) {

    const G4VEMDataSet* set = (zFluoDataVector[index])->GetComponent(j);
    G4double cross    = set->FindValue(kineticEnergy);

    p.push_back(cross);
    norm += cross;
  }

  if(norm == 0.0) return 0;

  G4double q = norm*G4UniformRand();

  for (G4int i=0; i<nElements; i++) {

    if(p[i] > q) {
       Z = (G4int)((*theElementVector)[i]->GetZ());
       break;
    }
    q -= p[i];
  }

  return Z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::ComputeDEDX(
                                 const G4ParticleDefinition* aParticle,
                                 const G4Material* material, 
				 G4double kineticEnergy) 
{  
  G4Proton* theProton = G4Proton::Proton();
  G4AntiProton* theAntiProton = G4AntiProton::AntiProton();
  G4double dedx = 0.0 ;
    
  G4double tscaled = kineticEnergy*proton_mass_c2/(aParticle->GetPDGMass()) ; 
  charge  = aParticle->GetPDGCharge() ;
  
  if(charge>0.0) {
    if(tscaled > protonHighEnergy) {
      dedx=G4EnergyLossTables::GetDEDX(theProton,tscaled,material) ;

    } else {
      dedx=ProtonParametrisedDEDX(material,tscaled) ;
    }

  } else {
    if(tscaled > antiProtonHighEnergy) {
      dedx=G4EnergyLossTables::GetDEDX(theAntiProton,tscaled,material); 
      
    } else {
      dedx=AntiProtonParametrisedDEDX(material,tscaled) ;
    }
  }    
  dedx *= theIonEffChargeModel->TheValue(aParticle, material, kineticEnergy) ;
  
  return dedx ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::BarkasTerm(const G4Material* material,
  				                  G4double kineticEnergy) const
//Function to compute the Barkas term for protons:
//
//Ref. Z_1^3 effect in the stopping power of matter for charged particles
//     J.C Ashley and R.H.Ritchie
//     Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397 
//
{
  static double FTable[47][2] = { 
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
  
  G4double BarkasTerm = 0.0;
  G4double AMaterial = 0.0;
  G4double ZMaterial = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int numberOfElements = material->GetNumberOfElements();
  
  for (G4int i = 0; i<numberOfElements; i++) {

    AMaterial = (*theElementVector)[i]->GetA()*mole/g;
    ZMaterial = (*theElementVector)[i]->GetZ();
    
    G4double X = 137.0 * 137.0 * beta2 / ZMaterial;
  
    // Variables to compute L_1
    G4double Eta0Chi = 0.8;
    G4double EtaChi = Eta0Chi * ( 1.0 + 6.02*pow( ZMaterial,-1.19 ) );
    G4double W = ( EtaChi * pow( ZMaterial,1.0/6.0 ) ) / sqrt(X); 
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
    
    BarkasTerm += FunctionOfW /( sqrt(ZMaterial * X) * X);
  }

  BarkasTerm *= twopi_mc2_rcl2 * (material->GetElectronDensity()) / beta2 ;

  return BarkasTerm;
}
        
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::BlochTerm(const G4Material* material,
                                                 G4double kineticEnergy,
                                                 G4double cSquare) const
//Function to compute the Bloch term for protons:
//
//Ref. Z_1^3 effect in the stopping power of matter for charged particles
//     J.C Ashley and R.H.Ritchie
//     Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397 
//
{
  G4double eloss = 0.0 ;
  G4double gamma = 1.0 + kineticEnergy / proton_mass_c2 ;
  G4double beta2 = 1.0 - 1.0/(gamma*gamma) ;
  G4double y = cSquare / (137.0*137.0*beta2) ;

  if(y < 0.05) {
    eloss = 1.202 ;
     
  } else {
    eloss = 1.0 / (1.0 + y) ;
    G4double de = eloss ;

    for(G4int i=2; de>eloss*0.01; i++) {
      de = 1.0/( i * (i*i + y)) ;
      eloss += de ;
    }
  }
  eloss *= -1.0 * y * cSquare * twopi_mc2_rcl2 * 
            (material->GetElectronDensity()) / beta2 ;
 
  return eloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4hLowEnergyIonisation::ElectronicLossFluctuation(
                                 const G4DynamicParticle* particle,
                                 const G4Material* material,
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

  G4int    imaterial   = material->GetIndex() ; 
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
  G4double threshold = G4std::max(deltaCutInKineticEnergyNow,ipotFluct);
  G4double loss, siga;

  G4double rmass = electron_mass_c2/particleMass;
  G4double tau   = tkin/particleMass; 
  G4double tau1 = tau+1.0; 
  G4double tau2 = tau*(tau+2.);
  G4double tmax    = 2.*electron_mass_c2*tau2/(1.+2.*tau1*rmass+rmass*rmass);

  
  if(tmax > threshold) tmax = threshold;
  G4double beta2 = tau2/(tau1*tau1);

  // Gaussian fluctuation 
  if(meanLoss > kappa*tmax || tmax < kappa*ipotFluct )
  {
    siga = tmax * (1.0-0.5*beta2) * step * twopi_mc2_rcl2 
         * electronDensity / beta2 ;

    // High velocity or negatively charged particle
    if( beta2 > 3.0*theBohrBeta2*zeff || charge < 0.0) {
      siga = sqrt( siga * chargeSquare ) ;

    // Low velocity - additional ion charge fluctuations according to
    // Q.Yang et al., NIM B61(1991)149-155.
    } else {
      G4double chu = theIonChuFluctuationModel->TheValue(particle, material);
      G4double yang = theIonYangFluctuationModel->TheValue(particle, material);
      siga = sqrt( siga * (chargeSquare * chu + yang)) ;
    }
    
    do {
        loss = G4RandGauss::shoot(meanLoss,siga);
    } while (loss < 0. || loss > 2.0*meanLoss);
    return loss;
  }

  // Non Gaussian fluctuation 
  static const G4double probLim = 0.01 ;
  static const G4double sumaLim = -log(probLim) ;
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

  w1 = tmax/ipotFluct;
  w2 = log(2.*electron_mass_c2*tau2);

  C = meanLoss*(1.-rateFluct)/(w2-ipotLogFluct-beta2);

  a1 = C*f1Fluct*(w2-e1LogFluct-beta2)/e1Fluct;
  a2 = C*f2Fluct*(w2-e2LogFluct-beta2)/e2Fluct;
  a3 = rateFluct*meanLoss*(tmax-ipotFluct)/(ipotFluct*tmax*log(w1));
  if(a1 < 0.0) a1 = 0.0;
  if(a2 < 0.0) a2 = 0.0;
  if(a3 < 0.0) a3 = 0.0;

  suma = a1+a2+a3;

  loss = 0.;


  if(suma < sumaLim)             // very small Step
    {
      e0 = material->GetIonisation()->GetEnergy0fluct();

      if(tmax == ipotFluct)
      {
        a3 = meanLoss/e0;

        if(a3>alim)
        {
          siga=sqrt(a3) ;
          p3 = G4std::max(0,G4int(G4RandGauss::shoot(a3,siga)+0.5));
        }
        else
          p3 = G4Poisson(a3);

        loss = p3*e0 ;

        if(p3 > 0)
          loss += (1.-2.*G4UniformRand())*e0 ;

      }
      else
      {
        tmax = tmax-ipotFluct+e0 ;
        a3 = meanLoss*(tmax-e0)/(tmax*e0*log(tmax/e0));

        if(a3>alim)
        {
          siga=sqrt(a3) ;
          p3 = G4std::max(0,int(G4RandGauss::shoot(a3,siga)+0.5));
        }
        else
          p3 = G4Poisson(a3);

        if(p3 > 0)
        {
          w = (tmax-e0)/tmax ;
          if(p3 > nmaxCont2)
          {
            dp3 = G4float(p3) ;
            corrfac = dp3/G4float(nmaxCont2) ;
            p3 = nmaxCont2 ;
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
        siga=sqrt(a1) ;
        p1 = G4std::max(0,G4int(G4RandGauss::shoot(a1,siga)+0.5));
      }
      else
       p1 = G4Poisson(a1);

      // excitation type 2
      if(a2>alim)
      {
        siga=sqrt(a2) ;
        p2 = G4std::max(0,G4int(G4RandGauss::shoot(a2,siga)+0.5));
      }
      else
        p2 = G4Poisson(a2);

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
        siga=sqrt(a3) ;
        p3 = G4std::max(0,G4int(G4RandGauss::shoot(a3,siga)+0.5));
      }
      else
        p3 = G4Poisson(a3);

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
            alfa1  = alfa*log(alfa)/(alfa-1.);
            ea     = na*ipotFluct*alfa1;
            sea    = ipotFluct*sqrt(na*(alfa-alfa1*alfa1));
            lossc += G4RandGauss::shoot(ea,sea);
          }
        }

        nb = G4int(G4float(p3)-na);
        if (nb > 0)
        {
          w2 = alfa*ipotFluct;
          w  = (tmax-w2)/tmax;      
          for (G4int k=0; k<nb; k++) lossc += w2/(1.-w*G4UniformRand());
        }
      }        
      loss += lossc;  
     }
    } 

  return loss ;
}
        
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::SetCutForSecondaryPhotons(G4double cut)
{
  minGammaEnergy = cut;
  deexcitationManager.SetCutForSecondaryPhotons(cut);
  theFluo = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::SetCutForAugerElectrons(G4double cut)
{
  minElectronEnergy = cut;
  deexcitationManager.SetCutForAugerElectrons(cut);
  theFluo = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::ActivateFluorescence(G4bool val)
{
  theFluo = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::ActivateAugerElectronProduction(G4bool val)
{
  deexcitationManager.ActivateAugerElectronProduction(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hLowEnergyIonisation::PrintInfoDefinition() const
{
  G4String comments = "  Knock-on electron cross sections . ";
  comments += "\n        Good description above the mean excitation energy.\n";
  comments += "        Delta ray energy sampled from  differential Xsection.";
  
  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << LowestKineticEnergy / eV << " eV " 
         << " to " << HighestKineticEnergy / TeV << " TeV "
         << " in " << TotBin << " bins."
 << "\n        Electronic stopping power model is  " 
 << theProtonTable
         << "\n        from " << protonLowEnergy / keV << " keV "
         << " to " << protonHighEnergy / MeV << " MeV " << "." << G4endl ;
  G4cout << "\n        Parametrisation model for antiprotons is  "
         << theAntiProtonTable
         << "\n        from " << antiProtonLowEnergy / keV << " keV "
         << " to " << antiProtonHighEnergy / MeV << " MeV " << "." << G4endl ;
  if(theBarkas){
  G4cout << "        Parametrization of the Barkas effect is switched on." 
         << G4endl ;
  }
  if(nStopping) {
  G4cout << "        Nuclear stopping power model is " << theNuclearTable 
         << G4endl ; 
  }

  G4bool printHead = true;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();
  
  // loop for materials 
  
  for (G4int j=0 ; j < numOfMaterials; j++) { 

    const G4Material* material= (*theMaterialTable)[j];
    G4double deltaCutNow = cutForDelta[(material->GetIndex())] ;   
    G4double eexc = material->GetIonisation()->GetMeanExcitationEnergy();

    if(eexc > deltaCutNow) {
      if(printHead) {
        printHead = false ;

        G4cout << "           material       min.delta energy(keV) " << G4endl;
        G4cout << G4endl;
      }

      G4cout << G4std::setw(20) << material->GetName()
	     << G4std::setw(15) << eexc/keV << G4endl;
    }
  }    
}
