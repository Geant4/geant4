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
// $Id: G4PenelopeCompton.cc,v 1.36 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// --------
// 12 Feb 2003   MG Pia       const argument in SelectRandomAtomForCompton
//                            Migration to "cuts per region"
// 14 Feb 2003   MG Pia       Corrected compilation errors and warnings
//                            from SUN
//                            Modified some variables to lowercase initial 
// 10 Mar 2003 V.Ivanchenko   Remove CutPerMaterial warning
// 13 Mar 2003 L.Pandola      Code "cleaned"  
// 20 Mar 2003 L.Pandola      ReadData() changed (performance improved) 
// 26 Mar 2003 L.Pandola      Added fluorescence
// 24 May 2003 MGP            Removed memory leak
// 09 Mar 2004 L.Pandola      Bug fixed in the generation of final state 
//                            (bug report # 585)
// 17 Mar 2004 L.Pandola      Removed unnecessary calls to std::pow(a,b)
// 18 Mar 2004 L.Pandola      Use of std::map (code review)
// 26 Mar 2008 L.Pandola      Add boolean flag to control atomic de-excitation
// 27 Mar 2008 L.Pandola      Re-named some variables to improve readability, 
//                            and check for strict energy conservation
// 03 Jun 2008 L.Pandola      Added further protection against non-conservation 
//                            of energy: it may happen because ionization energy
//                            from de-excitation manager and from Penelope internal 
//                            database do not match (difference is <10 eV, but may 
//                            give a e- with negative kinetic energy). 
//
// -------------------------------------------------------------------

#include "G4PenelopeCompton.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ForceCondition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4DynamicParticle.hh"
#include "G4VParticleChange.hh"
#include "G4ThreeVector.hh"
#include "G4EnergyLossTables.hh"
#include "G4VCrossSectionHandler.hh"
#include "G4CrossSectionHandler.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VRangeTest.hh"
#include "G4RangeTest.hh"
#include "G4ProductionCutsTable.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4AtomicDeexcitation.hh"
#include "G4PenelopeIntegrator.hh"
#include "G4MaterialCutsCouple.hh"


G4PenelopeCompton::G4PenelopeCompton(const G4String& processName)
  : G4VDiscreteProcess(processName),
    lowEnergyLimit(250*eV),
    highEnergyLimit(100*GeV),
    intrinsicLowEnergyLimit(10*eV),
    intrinsicHighEnergyLimit(100*GeV),
    energyForIntegration(0.0),
    ZForIntegration(1),
    nBins(200),
    cutForLowEnergySecondaryPhotons(250.0*eV),
    fUseAtomicDeexcitation(true)
{
  if (lowEnergyLimit < intrinsicLowEnergyLimit ||
      highEnergyLimit > intrinsicHighEnergyLimit)
    {
      G4Exception("G4PenelopeCompton::G4PenelopeCompton - energy outside intrinsic process validity range");
    }

  meanFreePathTable = 0;
  ionizationEnergy = new std::map<G4int,G4DataVector*>;
  hartreeFunction  = new std::map<G4int,G4DataVector*>;
  occupationNumber = new std::map<G4int,G4DataVector*>;

  rangeTest = new G4RangeTest;

  ReadData(); //Read data from file

  if (verboseLevel > 0)
    {
      G4cout << GetProcessName() << " is created " << G4endl
	     << "Energy range: "
	     << lowEnergyLimit / keV << " keV - "
	     << highEnergyLimit / GeV << " GeV"
	     << G4endl;
    }


   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4PenelopeCompton is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;

}

G4PenelopeCompton::~G4PenelopeCompton()
{
  delete meanFreePathTable;
  delete rangeTest;

  for (size_t i=0;i<matCrossSections->size();i++)
    {
      delete (*matCrossSections)[i];
    }

  delete matCrossSections;

  for (G4int Z=1;Z<100;Z++)
    {
      if (ionizationEnergy->count(Z)) delete (ionizationEnergy->find(Z)->second);
      if (hartreeFunction->count(Z)) delete (hartreeFunction->find(Z)->second);
      if (occupationNumber->count(Z)) delete (occupationNumber->find(Z)->second);
    }
  delete ionizationEnergy;
  delete hartreeFunction;
  delete occupationNumber;
}

void G4PenelopeCompton::BuildPhysicsTable(const G4ParticleDefinition& )
{
  G4DataVector energyVector;
  G4double dBin = std::log10(highEnergyLimit/lowEnergyLimit)/nBins;
  for (G4int i=0;i<nBins;i++)
    {
      energyVector.push_back(std::pow(10.,std::log10(lowEnergyLimit)+i*dBin));
    }

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  G4int nMaterials = G4Material::GetNumberOfMaterials();
  G4VDataSetAlgorithm* algo = new G4LogLogInterpolation();

  //size_t nOfBins = energyVector.size();
  //size_t bin=0;

  G4DataVector* energies;
  G4DataVector* data;

  matCrossSections = new std::vector<G4VEMDataSet*>;

  for (G4int m=0; m<nMaterials; m++)
    {
      const G4Material* material= (*materialTable)[m];
      G4int nElements = material->GetNumberOfElements();
      const G4ElementVector* elementVector = material->GetElementVector();
      const G4double* nAtomsPerVolume = material->GetAtomicNumDensityVector();

      G4VEMDataSet* setForMat = new G4CompositeEMDataSet(algo,1.,1.);

      for (G4int i=0; i<nElements; i++) {
 
        G4int Z = (G4int) (*elementVector)[i]->GetZ();
        G4double density = nAtomsPerVolume[i];
	G4double cross=0.0;
        energies = new G4DataVector;
        data = new G4DataVector;


        for (size_t bin=0; bin<energyVector.size(); bin++)
	  {
	    G4double e = energyVector[bin];
	    energies->push_back(e);
	    cross = density * CrossSection(e,Z); 
	    data->push_back(cross);
	  }

        G4VEMDataSet* elSet = new G4EMDataSet(i,energies,data,algo,1.,1.);
        setForMat->AddComponent(elSet);
      }

      matCrossSections->push_back(setForMat);
    }


  //Build the mean free path table! 
  G4double matCS = 0.0;
  G4VEMDataSet* matCrossSet = new G4CompositeEMDataSet(algo,1.,1.);
  G4VEMDataSet* materialSet = new G4CompositeEMDataSet(algo,1.,1.);
 
 
  for (G4int m=0; m<nMaterials; m++)
    { 
      energies = new G4DataVector;
      data = new G4DataVector;
      const G4Material* material= (*materialTable)[m];
      material= (*materialTable)[m];
      for (size_t bin=0; bin<energyVector.size(); bin++)
	{
	  G4double energy = energyVector[bin];
	  energies->push_back(energy);
	  matCrossSet = (*matCrossSections)[m]; 
	  matCS = 0.0;
          G4int nElm = matCrossSet->NumberOfComponents();
          for(G4int j=0; j<nElm; j++) {
            matCS += matCrossSet->GetComponent(j)->FindValue(energy);
	  }
	  if (matCS > 0.)
	    {
	      data->push_back(1./matCS);
	    }
	  else
	    {
	      data->push_back(DBL_MAX);
	    }
	}
      G4VEMDataSet* dataSet = new G4EMDataSet(m,energies,data,algo,1.,1.);   
      materialSet->AddComponent(dataSet);
    }
  meanFreePathTable = materialSet; 
}

G4VParticleChange* G4PenelopeCompton::PostStepDoIt(const G4Track& aTrack, 
						   const G4Step&  aStep)
{
  //Penelope model

  aParticleChange.Initialize(aTrack);
  
  // Dynamic particle quantities  
  const G4DynamicParticle* incidentPhoton = aTrack.GetDynamicParticle();
  G4double photonEnergy0 = incidentPhoton->GetKineticEnergy();

  if (photonEnergy0 <= lowEnergyLimit)
    {
      aParticleChange.ProposeTrackStatus(fStopAndKill);
      aParticleChange.ProposeEnergy(0.);
      aParticleChange.ProposeLocalEnergyDeposit(photonEnergy0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }

  G4ParticleMomentum photonDirection0 = incidentPhoton->GetMomentumDirection();

  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
  const G4Material* material = couple->GetMaterial();
  
  G4int Z = SelectRandomAtomForCompton(material,photonEnergy0);
  const G4int nmax = 64;
  G4double rn[nmax],pac[nmax];
  
  G4double ki,ki1,ki2,ki3,taumin,a1,a2;
  G4double tau,TST;
  G4double S=0.0;
  G4double epsilon,cosTheta;
  G4double harFunc = 0.0;
  G4int occupNb= 0;
  G4double ionEnergy=0.0;
  G4int nosc = occupationNumber->find(Z)->second->size();
  G4int iosc = nosc;
  ki = photonEnergy0/electron_mass_c2;
  ki2 = 2*ki+1.0;
  ki3 = ki*ki;
  ki1 = ki3-ki2-1.0;
  taumin = 1.0/ki2;
  a1 = std::log(ki2);
  a2 = a1+2.0*ki*(1.0+ki)/(ki2*ki2);
  //If the incoming photon is above 5 MeV, the quicker approach based on the 
  //pure Klein-Nishina formula is used
  if (photonEnergy0 > 5*MeV)
    {
      do{
	do{
	  if ((a2*G4UniformRand()) < a1)
	    {
	      tau = std::pow(taumin,G4UniformRand());
	    }
	  else
	    {
	      tau = std::sqrt(1.0+G4UniformRand()*(taumin*taumin-1.0));
	    }
	  //rejection function
	  TST = (1+tau*(ki1+tau*(ki2+tau*ki3)))/(ki3*tau*(1.0+tau*tau));
	}while (G4UniformRand()> TST);
	epsilon=tau;
	cosTheta = 1.0 - (1.0-tau)/(ki*tau);
	//Target shell electrons
	TST = Z*G4UniformRand();
	iosc = nosc;
	S=0.0;
	for (G4int j=0;j<nosc;j++)
	  {
	    occupNb = (G4int) (*(occupationNumber->find(Z)->second))[j];
	    S = S + occupNb;
	    if (S > TST) iosc = j;
	    if (S > TST) break; 
	  }
	ionEnergy = (*(ionizationEnergy->find(Z)->second))[iosc];
      }while((epsilon*photonEnergy0-photonEnergy0+ionEnergy) >0);
    }
  else //photonEnergy0<5 MeV
    {
      //Incoherent scattering function for theta=PI
      G4double s0=0.0;
      G4double pzomc=0.0,rni=0.0;
      G4double aux=0.0;
      for (G4int i=0;i<nosc;i++){
	ionEnergy = (*(ionizationEnergy->find(Z)->second))[i];
	if (photonEnergy0 > ionEnergy)
	  {
	    G4double aux = photonEnergy0*(photonEnergy0-ionEnergy)*2.0;
	    harFunc = (*(hartreeFunction->find(Z)->second))[i]/fine_structure_const;
	    occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
	    pzomc = harFunc*(aux-electron_mass_c2*ionEnergy)/
	       (electron_mass_c2*std::sqrt(2.0*aux+ionEnergy*ionEnergy));
	    if (pzomc > 0) 
	      {
		rni = 1.0-0.5*std::exp(0.5-(std::sqrt(0.5)+std::sqrt(2.0)*pzomc)*(std::sqrt(0.5)+std::sqrt(2.0)*pzomc));
	      }
	    else
	      {
		rni = 0.5*std::exp(0.5-(std::sqrt(0.5)-std::sqrt(2.0)*pzomc)*(std::sqrt(0.5)-std::sqrt(2.0)*pzomc));
	      }
	    s0 = s0 + occupNb*rni;
	  }
      }
      
      //Sampling tau
      G4double cdt1;
      do
	{
	  if ((G4UniformRand()*a2) < a1)
	    {
	      tau = std::pow(taumin,G4UniformRand());
	    }
	  else
	    {
	      tau = std::sqrt(1.0+G4UniformRand()*(taumin*taumin-1.0));
	    }
	  cdt1 = (1.0-tau)/(ki*tau);
	  S=0.0;
	  //Incoherent scattering function
	  for (G4int i=0;i<nosc;i++){
	    ionEnergy = (*(ionizationEnergy->find(Z)->second))[i];
	    if (photonEnergy0 > ionEnergy) //sum only on excitable levels
	      {
		aux = photonEnergy0*(photonEnergy0-ionEnergy)*cdt1;
		harFunc = (*(hartreeFunction->find(Z)->second))[i]/fine_structure_const;
		occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
		pzomc = harFunc*(aux-electron_mass_c2*ionEnergy)/
		  (electron_mass_c2*std::sqrt(2.0*aux+ionEnergy*ionEnergy));
		if (pzomc > 0) 
		  {
		    rn[i] = 1.0-0.5*std::exp(0.5-(std::sqrt(0.5)+std::sqrt(2.0)*pzomc)*
					     (std::sqrt(0.5)+std::sqrt(2.0)*pzomc));
		  }
		else
		  {
		    rn[i] = 0.5*std::exp(0.5-(std::sqrt(0.5)-std::sqrt(2.0)*pzomc)*
					 (std::sqrt(0.5)-std::sqrt(2.0)*pzomc));
		  }
		S = S + occupNb*rn[i];
		pac[i] = S;
	      }
	    else
	      {
		pac[i] = S-(1e-06);
	      }
	  }
	  //Rejection function
	  TST = S*(1.0+tau*(ki1+tau*(ki2+tau*ki3)))/(ki3*tau*(1.0+tau*tau));  
	}while ((G4UniformRand()*s0) > TST);

      //Target electron shell
      cosTheta = 1.0 - cdt1;
      G4double fpzmax=0.0,fpz=0.0;
      G4double A=0.0;
      do
	{
	  do
	    {
	      TST =S*G4UniformRand();
	      iosc=nosc;
	      for (G4int i=0;i<nosc;i++){
		if (pac[i]>TST) iosc = i;
		if (pac[i]>TST) break; 
	      }
	      A = G4UniformRand()*rn[iosc];
	      harFunc = (*(hartreeFunction->find(Z)->second))[iosc]/fine_structure_const;
	      occupNb = (G4int) (*(occupationNumber->find(Z)->second))[iosc];
	      if (A < 0.5) {
		pzomc = (std::sqrt(0.5)-std::sqrt(0.5-std::log(2.0*A)))/
		  (std::sqrt(2.0)*harFunc);
	      }
	      else
		{
		  pzomc = (std::sqrt(0.5-std::log(2.0-2.0*A))-std::sqrt(0.5))/
		    (std::sqrt(2.0)*harFunc);
		}
	    } while (pzomc < -1);
	  // F(EP) rejection
	  G4double XQC = 1.0+tau*(tau-2.0*cosTheta);
	  G4double AF = std::sqrt(XQC)*(1.0+tau*(tau-cosTheta)/XQC);
	  if (AF > 0) {
	    fpzmax = 1.0+AF*0.2;
	  }
	  else
	    {
	      fpzmax = 1.0-AF*0.2;
	    }
	  fpz = 1.0+AF*std::max(std::min(pzomc,0.2),-0.2);
	}while ((fpzmax*G4UniformRand())>fpz);
  
      //Energy of the scattered photon
      G4double T = pzomc*pzomc;
      G4double b1 = 1.0-T*tau*tau;
      G4double b2 = 1.0-T*tau*cosTheta;
      if (pzomc > 0.0)
	{
	  epsilon = (tau/b1)*(b2+std::sqrt(std::abs(b2*b2-b1*(1.0-T))));
	}
      else
	{
	  epsilon = (tau/b1)*(b2-std::sqrt(std::abs(b2*b2-b1*(1.0-T))));
	}
    }
  
  G4double sinTheta = std::sqrt(1-cosTheta*cosTheta);
  G4double phi = twopi * G4UniformRand() ;
  G4double dirx = sinTheta * std::cos(phi);
  G4double diry = sinTheta * std::sin(phi);
  G4double dirz = cosTheta ;

  // Update G4VParticleChange for the scattered photon 
  
  G4ThreeVector photonDirection1(dirx,diry,dirz);
  photonDirection1.rotateUz(photonDirection0);
  aParticleChange.ProposeMomentumDirection(photonDirection1) ;
  G4double photonEnergy1 = epsilon * photonEnergy0;   

  if (photonEnergy1 > 0.)
    {
      aParticleChange.ProposeEnergy(photonEnergy1) ;
    }
  else
    {    
      aParticleChange.ProposeEnergy(0.) ;
      aParticleChange.ProposeTrackStatus(fStopAndKill);
    }


  // Kinematics of the scattered electron   
  G4double diffEnergy = photonEnergy0*(1-epsilon);
  ionEnergy = (*(ionizationEnergy->find(Z)->second))[iosc];
  G4double Q2 = photonEnergy0*photonEnergy0+photonEnergy1*(photonEnergy1-2.0*photonEnergy0*cosTheta);
  G4double cosThetaE; //scattering angle for the electron
  if (Q2 > 1.0e-12)
    {
      cosThetaE = (photonEnergy0-photonEnergy1*cosTheta)/std::sqrt(Q2);
    }
  else
    {
      cosThetaE = 1.0;
    }
  G4double sinThetaE = std::sqrt(1-cosThetaE*cosThetaE);
  //initialize here, then check photons created by Atomic-Deexcitation, and the final state e-
  G4int nbOfSecondaries = 0; 
  
  std::vector<G4DynamicParticle*>* photonVector=0;

  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  const G4AtomicShell* shell = transitionManager->Shell(Z,iosc);
  G4double bindingEnergy = shell->BindingEnergy();
  G4int shellId = shell->ShellId();
  G4double ionEnergyInPenelopeDatabase = ionEnergy;
  ionEnergy = std::max(bindingEnergy,ionEnergyInPenelopeDatabase); //protection against energy non-conservation 

  G4double eKineticEnergy = diffEnergy - ionEnergy; //subtract the excitation energy. If not emitted by fluorescence,
  //the ionization energy is deposited as local energy deposition
  G4double localEnergyDeposit = ionEnergy; 
  G4double energyInFluorescence = 0.; //testing purposes only

  if (eKineticEnergy < 0) 
    {
      //It means that there was some problem/mismatch between the two databases. Try to make it work
      //In this case available Energy (diffEnergy) < ionEnergy
      //Full residual energy is deposited locally
      localEnergyDeposit = diffEnergy;
      eKineticEnergy = 0.0;
    }

  //the local energy deposit is what remains: part of this may be spent for fluorescence.
 
  if (fUseAtomicDeexcitation)
    { 
      G4int nPhotons=0;
      
      const G4ProductionCutsTable* theCoupleTable=
	G4ProductionCutsTable::GetProductionCutsTable();
      size_t indx = couple->GetIndex();

      G4double cutg = (*(theCoupleTable->GetEnergyCutsVector(0)))[indx];
      cutg = std::max(cutForLowEnergySecondaryPhotons,cutg);

      G4double cute = (*(theCoupleTable->GetEnergyCutsVector(1)))[indx];
      cute = std::max(cutForLowEnergySecondaryPhotons,cute);
     
      G4DynamicParticle* aPhoton;
      G4AtomicDeexcitation deexcitationManager;
      
      if (Z>5 && (localEnergyDeposit > cutg || localEnergyDeposit > cute))
	{
	  photonVector = deexcitationManager.GenerateParticles(Z,shellId);
	  for (size_t k=0;k<photonVector->size();k++){
	    aPhoton = (*photonVector)[k];
	    if (aPhoton)
	      {
		G4double itsCut = cutg;
		if (aPhoton->GetDefinition() == G4Electron::Electron()) itsCut = cute;
		G4double itsEnergy = aPhoton->GetKineticEnergy();
		if (itsEnergy > itsCut && itsEnergy <= localEnergyDeposit)
		  {
		    nPhotons++;
		    localEnergyDeposit -= itsEnergy;
		    energyInFluorescence += itsEnergy;
		  }
		else
		  {
		    delete aPhoton;
		    (*photonVector)[k]=0;
		  }
	      }
	  }
	}
      nbOfSecondaries=nPhotons;
    }

  
  // Generate the electron only if with large enough range w.r.t. cuts and safety 
  G4double safety = aStep.GetPostStepPoint()->GetSafety();
  G4DynamicParticle* electron = 0;
  if (rangeTest->Escape(G4Electron::Electron(),couple,eKineticEnergy,safety) && 
      eKineticEnergy>cutForLowEnergySecondaryPhotons)
    {
      G4double xEl = sinThetaE * std::cos(phi+pi); 
      G4double yEl = sinThetaE * std::sin(phi+pi);
      G4double zEl = cosThetaE;
      G4ThreeVector eDirection(xEl,yEl,zEl); //electron direction
      eDirection.rotateUz(photonDirection0);
      electron = new G4DynamicParticle (G4Electron::Electron(),
					eDirection,eKineticEnergy) ;
      nbOfSecondaries++;
    }
  else
    {
      localEnergyDeposit += eKineticEnergy;
    }

  aParticleChange.SetNumberOfSecondaries(nbOfSecondaries);
  if (electron) aParticleChange.AddSecondary(electron);

  //This block below is executed only if there is at least one secondary photon produced by
  //AtomicDeexcitation
  if (photonVector)
    {
      for (size_t ll=0;ll<photonVector->size();ll++)
	{
	  if ((*photonVector)[ll]) aParticleChange.AddSecondary((*photonVector)[ll]);
	}
    }
  delete photonVector;
  if (localEnergyDeposit < 0)
    {
      G4cout << "WARNING-" 
	     << "G4PenelopeCompton::PostStepDoIt - Negative energy deposit"
	     << G4endl;
      localEnergyDeposit=0.;
    }
  aParticleChange.ProposeLocalEnergyDeposit(localEnergyDeposit);
  

  if (verboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4PenelopeCompton" << G4endl;
      G4cout << "Incoming photon energy: " << photonEnergy0/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Scattered photon: " << photonEnergy1/keV << " keV" << G4endl;
      G4double electronEnergy = 0.;
      if (electron)
	electronEnergy = eKineticEnergy;
      G4cout << "Scattered electron " << electronEnergy/keV << " keV" << G4endl;
      G4cout << "Fluorescence: " << energyInFluorescence/keV << " keV" << G4endl;
      G4cout << "Local energy deposit " << localEnergyDeposit/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (photonEnergy1+electronEnergy+energyInFluorescence+localEnergyDeposit)/keV << 
	" keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }

  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
}

G4bool G4PenelopeCompton::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() ); 
}

G4double G4PenelopeCompton::GetMeanFreePath(const G4Track& track, 
					    G4double, // previousStepSize
					    G4ForceCondition*)
{
  const G4DynamicParticle* photon = track.GetDynamicParticle();
  G4double energy = photon->GetKineticEnergy();
  G4Material* material = track.GetMaterial();
  size_t materialIndex = material->GetIndex();

  G4double meanFreePath;
  if (energy > highEnergyLimit) meanFreePath = meanFreePathTable->FindValue(highEnergyLimit,materialIndex);
  else if (energy < lowEnergyLimit) meanFreePath = DBL_MAX;
  else meanFreePath = meanFreePathTable->FindValue(energy,materialIndex);
  return meanFreePath;
}


void G4PenelopeCompton::ReadData()
{
  char* path = getenv("G4LEDATA");
  if (!path)
    {
      G4String excep = "G4PenelopeCompton - G4LEDATA environment variable not set!";
      G4Exception(excep);
    }
  G4String pathString(path);
  G4String pathFile = pathString + "/penelope/compton-pen.dat";
  std::ifstream file(pathFile);
  std::filebuf* lsdp = file.rdbuf();
  
  if (!(lsdp->is_open()))
    {
      G4String excep = "G4PenelopeCompton - data file " + pathFile + " not found!";
      G4Exception(excep);
    }

  G4int k1,test,test1;
  G4double a1,a2;
  G4int Z=1,nLevels=0;
  G4DataVector* f;
  G4DataVector* u;
  G4DataVector* j;

  do{
    f = new G4DataVector;
    u = new G4DataVector;
    j = new G4DataVector;
    file >> Z >> nLevels;
    for (G4int h=0;h<nLevels;h++){
      file >> k1 >> a1 >> a2;
      f->push_back((G4double) k1);
      u->push_back(a1);
      j->push_back(a2);
    }
    ionizationEnergy->insert(std::make_pair(Z,u));
    hartreeFunction->insert(std::make_pair(Z,j));
    occupationNumber->insert(std::make_pair(Z,f));
    file >> test >> test1; //-1 -1 close the data for each Z
    if (test > 0) {
      G4String excep = "G4PenelopeCompton - data file corrupted!";
      G4Exception(excep);
    }
  }while (test != -2); //the very last Z is closed with -2 instead of -1
}

G4double G4PenelopeCompton::CrossSection(G4double energy,G4int Z)
{
  G4double cs=0.0;
  energyForIntegration=energy; 
  ZForIntegration = Z;
  if (energy< 5*MeV)
    {
      G4PenelopeIntegrator<G4PenelopeCompton,G4double (G4PenelopeCompton::*)(G4double)> theIntegrator;
      cs = theIntegrator.Calculate(this,&G4PenelopeCompton::DifferentialCrossSection,-1.0,1.0,1e-05);
    }
  else
    {
      G4double ki=energy/electron_mass_c2;
      G4double ki3=ki*ki;
      G4double ki2=1.0+2*ki;
      G4double ki1=ki3-ki2-1.0;
      G4double t0=1.0/(ki2);
      G4double csl = 0.5*ki3*t0*t0+ki2*t0+ki1*std::log(t0)-(1.0/t0);
      G4int nosc = occupationNumber->find(Z)->second->size();
      for (G4int i=0;i<nosc;i++)
	{
	  G4double ionEnergy = (*(ionizationEnergy->find(Z)->second))[i];
	  G4double tau=(energy-ionEnergy)/energy;
	  if (tau > t0)
	    {
	      G4double csu = 0.5*ki3*tau*tau+ki2*tau+ki1*std::log(tau)-(1.0/tau);
	      G4int f = (G4int) (*(occupationNumber->find(Z)->second))[i];
	      cs = cs + f*(csu-csl);
	    }
	}
      cs=pi*classic_electr_radius*classic_electr_radius*cs/(ki*ki3);
    }
  return cs;
}

  
G4double G4PenelopeCompton::DifferentialCrossSection(G4double cosTheta)
{
  const G4double k2 = std::sqrt(2.0);
  const G4double k1 = std::sqrt(0.5);
  const G4double k12 = 0.5;
  G4double cdt1 = 1.0-cosTheta;
  G4double energy = energyForIntegration;
  G4int Z = ZForIntegration;
  G4double ionEnergy=0.0,Pzimax=0.0,XKN=0.0;
  G4double diffCS=0.0;
  G4double x=0.0,siap=0.0;
  G4double harFunc=0.0;
  G4int occupNb;
  //energy of Compton line;
  G4double EOEC = 1.0+(energy/electron_mass_c2)*cdt1; 
  G4double ECOE = 1.0/EOEC;
  //Incoherent scattering function (analytical profile)
  G4double sia = 0.0;
  G4int nosc = occupationNumber->find(Z)->second->size();
  for (G4int i=0;i<nosc;i++){
    ionEnergy = (*(ionizationEnergy->find(Z)->second))[i];
    //Sum only of those shells for which E>Eion
    if (energy > ionEnergy)
      {
	G4double aux = energy * (energy-ionEnergy)*cdt1;
	Pzimax = (aux - electron_mass_c2*ionEnergy)/(electron_mass_c2*std::sqrt(2*aux+ionEnergy*ionEnergy));
	harFunc = (*(hartreeFunction->find(Z)->second))[i]/fine_structure_const;
	occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
	x = harFunc*Pzimax;
	if (x > 0) 
	  {
	    siap = 1.0-0.5*std::exp(k12-(k1+k2*x)*(k1+k2*x));
	  }
	else
	  {
	    siap = 0.5*std::exp(k12-(k1-k2*x)*(k1-k2*x));
	  }
	sia = sia + occupNb*siap; //sum of all contributions;
      }
  }
  XKN = EOEC+ECOE-1+cosTheta*cosTheta;
  diffCS = pi*classic_electr_radius*classic_electr_radius*ECOE*ECOE*XKN*sia;
  return diffCS;
}

G4int G4PenelopeCompton::SelectRandomAtomForCompton(const G4Material* material,G4double energy) const
{
  G4int nElements = material->GetNumberOfElements();
  //Special case: the material consists of one element
  if (nElements == 1)
    {
      G4int Z = (G4int) material->GetZ();
      return Z;
    }

  //Composite material
  const G4ElementVector* elementVector = material->GetElementVector();
  size_t materialIndex = material->GetIndex();

  G4VEMDataSet* materialSet = (*matCrossSections)[materialIndex]; 
  G4double materialCrossSection0 = 0.0;
  G4DataVector cross;
  cross.clear();
  G4int i;
  for (i=0;i<nElements;i++)
    {
      G4double cr = (materialSet->GetComponent(i))->FindValue(energy);
      materialCrossSection0 += cr;
      cross.push_back(materialCrossSection0); //cumulative cross section
    }

  G4double random = G4UniformRand()*materialCrossSection0;
  for (i=0;i<nElements;i++)
    {
      if (random <= cross[i]) return (G4int) (*elementVector)[i]->GetZ();
    }
  //It should never get here
  return 0;
}
  
