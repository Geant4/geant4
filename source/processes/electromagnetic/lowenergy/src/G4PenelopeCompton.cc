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
// $Id: G4PenelopeCompton.cc,v 1.3 2002-12-12 09:55:33 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// --------
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
#include "G4PenelopeIntegrator.hh"

#include "G4CutsPerMaterialWarning.hh"

G4PenelopeCompton::G4PenelopeCompton(const G4String& processName)
  : G4VDiscreteProcess(processName),
    lowEnergyLimit(250*eV),             
    highEnergyLimit(100*GeV),
    intrinsicLowEnergyLimit(10*eV),
    intrinsicHighEnergyLimit(100*GeV),
    energyForIntegration(0.0),
    ZForIntegration(1),
    nBins(200)
{
  if (lowEnergyLimit < intrinsicLowEnergyLimit || 
      highEnergyLimit > intrinsicHighEnergyLimit)
    {
      G4Exception("G4PenelopeCompton::G4PenelopeCompton - energy outside intrinsic process validity range");
    }

  meanFreePathTable = 0;
  IonizationEnergy = new G4std::vector<G4DataVector*>;
  HartreeFunction  = new G4std::vector<G4DataVector*>;
  OccupationNumber = new G4std::vector<G4DataVector*>;
  rangeTest = new G4RangeTest;
 
  ReadData(); //Legge i dati da file!
 
  if (verboseLevel > 0) 
    {
      G4cout << GetProcessName() << " is created " << G4endl
	     << "Energy range: " 
	     << lowEnergyLimit / keV << " keV - "
	     << highEnergyLimit / GeV << " GeV" 
	     << G4endl;
    }
}
 
G4PenelopeCompton::~G4PenelopeCompton()
{
  delete meanFreePathTable;
  delete rangeTest;
  delete matCrossSections;
  delete IonizationEnergy;
  delete HartreeFunction;
  delete OccupationNumber;
}

void G4PenelopeCompton::BuildPhysicsTable(const G4ParticleDefinition& photon)
{
  
  G4CutsPerMaterialWarning warning;
  warning.PrintWarning(&photon);

  G4DataVector energyVector;
  G4double dBin = log10(highEnergyLimit/lowEnergyLimit)/nBins;
  for (G4int i=0;i<nBins;i++)
    {
      energyVector.push_back(pow(10.,log10(lowEnergyLimit)+i*dBin));
    }

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  G4int nMaterials = G4Material::GetNumberOfMaterials();
  G4VDataSetAlgorithm* algo = new G4LogLogInterpolation();

  size_t nOfBins = energyVector.size();
  size_t bin=0;

  G4DataVector* energies;
  G4DataVector* data;

  matCrossSections = new G4std::vector<G4VEMDataSet*>;

  

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


        for (bin=0; bin<nOfBins; bin++)
	  {
	    G4double e = energyVector[bin];
	    energies->push_back(e);
	    cross = density * CrossSection(e,Z); //sezione d'urto
	    data->push_back(cross);
	  }

        G4VEMDataSet* elSet = new G4EMDataSet(i,energies,data,algo,1.,1.);
        setForMat->AddComponent(elSet);
      }

      matCrossSections->push_back(setForMat);
    }

 //sono state calcolate le sezioni d'urto dei vari materiali
 G4double matCS = 0.0;
 G4VEMDataSet* matCrossSet = new G4CompositeEMDataSet(algo,1.,1.);
 G4VEMDataSet* materialSet = new G4CompositeEMDataSet(algo,1.,1.);
 
 for (m=0; m<nMaterials; m++)
   { 
      energies = new G4DataVector;
      data = new G4DataVector;
      material= (*materialTable)[m];
      for (bin=0; bin<nOfBins; bin++)
	{
	  G4double energy = energyVector[bin];
	  energies->push_back(energy);
	  matCrossSet = (*matCrossSections)[m]; //E'un composite!
          G4int nElm = matCrossSet->NumberOfComponents();//una per ogni elemento
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
      aParticleChange.SetStatusChange(fStopAndKill);
      aParticleChange.SetEnergyChange(0.);
      aParticleChange.SetLocalEnergyDeposit(photonEnergy0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }

  G4ParticleMomentum photonDirection0 = incidentPhoton->GetMomentumDirection();

  G4Material* material = aTrack.GetMaterial();
  G4int Z = SelectRandomAtomForCompton(material,photonEnergy0);
  //La routine di estrazione funziona
  const G4int nmax = 64;
  G4double rn[nmax],pac[nmax];
  
  G4double ki,ki1,ki2,ki3,taumin,a1,a2;
  G4double tau,TST;
  G4double S=0.0;
  G4double epsilon,cosTheta;
  G4double HarFunc = 0.0;
  G4int OccupNb= 0;
  G4double IonEnergy=0.0;
  G4int nosc = ((*OccupationNumber)[Z-1])->size();
  G4int iosc = nosc;
  ki = photonEnergy0/electron_mass_c2;
  ki2 = 2*ki+1.0;
  ki3 = ki*ki;
  ki1 = ki3-ki2-1.0;
  taumin = 1.0/ki2;
  a1 = log(ki2);
  a2 = a1+2.0*ki*(1.0+ki)/(ki2*ki2);
  if (photonEnergy0 > 5*MeV)
    {
      do{
	do{
	  if ((a2*G4UniformRand()) < a1)
	    {
	      tau = taumin*G4UniformRand();
	    }
	  else
	    {
	      tau = sqrt(1.0+G4UniformRand()*(taumin*taumin-1.0));
	    }
	  //rejection function
	  TST = (1+tau*(ki1+tau*(ki2+tau*ki3)))/(ki3*tau*(1.0+tau*tau));
	}while (G4UniformRand()> TST);
	epsilon=tau;
	cosTheta = 1.0 - (1.0-tau)/(ki*tau);
	//Target shell electron
	TST = Z*G4UniformRand(); 
	for (G4int j=0;j<nosc;j++)
	  {
	    OccupNb = (G4int) (*((*OccupationNumber)[Z-1]))[j];
	    S = S + OccupNb;
	    if (S > TST) iosc = j;
	    if (S > TST) break; //funziona anche dentro l'if??
	  }
	IonEnergy = (*((*IonizationEnergy)[Z-1]))[iosc];
      }while((epsilon*photonEnergy0-photonEnergy0+IonEnergy) >0);
    }

  else //photonEnergy0<5 MeV
    {
      //Incoherent scattering function for theta=PI
      G4double s0=0.0;
      G4double pzomc=0.0,rni=0.0;
      G4double aux=0.0;
      for (G4int i=0;i<nosc;i++){
	IonEnergy = (*((*IonizationEnergy)[Z-1]))[i];
	if (IonEnergy < photonEnergy0)
	  {
	    G4double aux = photonEnergy0*(photonEnergy0-IonEnergy)*2.0;
	    HarFunc = (*((*HartreeFunction)[Z-1]))[i];
	    OccupNb = (G4int) (*((*OccupationNumber)[Z-1]))[i];
	    pzomc = HarFunc*(aux-electron_mass_c2*IonEnergy)/(electron_mass_c2*sqrt(2.0*aux+pow(IonEnergy,2)));
	    if (pzomc > 0) 
	      {
		rni = 1.0-0.5*exp(0.5-pow(sqrt(0.5)+sqrt(2.0)*pzomc,2));
	      }
	    else
	      {
		rni = 0.5*exp(0.5-pow(sqrt(0.5)-sqrt(2.0)*pzomc,2));
	      }
	    s0 = s0 + OccupNb*rni;
	  }
      }
      
      //Sampling tau
      G4double cdt1;
      do
	{
	  if ((G4UniformRand()*a2) > a1)
	    {
	      tau = pow(taumin,G4UniformRand());
	    }
	  else
	    {
	      tau = sqrt(1.0+G4UniformRand()*(taumin*taumin-1.0));
	    }
	  cdt1 = (1.0-tau)/(ki*tau);
	  
	  //Incoherent scattering function
	  for (i=0;i<nosc;i++){
	    IonEnergy = (*((*IonizationEnergy)[Z-1]))[i];
	    if (IonEnergy < photonEnergy0)
	      {
		aux = photonEnergy0*(photonEnergy0-IonEnergy)*cdt1;
		HarFunc = (*((*HartreeFunction)[Z-1]))[i];
		OccupNb = (G4int) (*((*OccupationNumber)[Z-1]))[i];
		pzomc = HarFunc*(aux-electron_mass_c2*IonEnergy)/(electron_mass_c2*sqrt(2.0*aux+pow(IonEnergy,2)));
		if (pzomc > 0) 
		  {
		    rn[i] = 1.0-0.5*exp(0.5-pow(sqrt(0.5)+sqrt(2.0)*pzomc,2));
		  }
		else
		  {
		    rn[i] = 0.5*exp(0.5-pow(sqrt(0.5)-sqrt(2.0)*pzomc,2));
		  }
		S = S + OccupNb*rn[i];
		pac[i] = S;
	      }
	    else
	      {
		pac[i] = S-1e-06;
	      }
	  }
	  //Rejection function
	  TST = S*(1.0+tau*(ki1+tau*(ki2+tau*ki3)))/(ki3*tau*(1.0+tau*tau));  
	}while ((G4UniformRand()*s0) > TST);
      
      //Target electron shell
      cosTheta = 1.0 - cdt1;
      iosc=nosc;
      G4double fpzmax=0.0,fpz=0.0;
      do
	{
	  do
	    {
	      TST =S*G4UniformRand();
	      for (i=0;i<nosc;i++){
		if (pac[i]>TST) iosc = i;
		if (pac[i]>TST) break; //funziona anche nello stesso if?
	      }
	      G4double A = G4UniformRand()*rn[iosc];
	      HarFunc = (*((*HartreeFunction)[Z-1]))[iosc];
	      OccupNb = (G4int) (*((*OccupationNumber)[Z-1]))[iosc];
	      if (A < 0.5) {
		pzomc = (sqrt(0.5)-sqrt(0.5-log(2*A)))/(sqrt(2)*HarFunc);
	      }
	      else
		{
		  pzomc = (sqrt(0.5-log(2.0-2.0*A))-sqrt(0.5))/(sqrt(2.0)*HarFunc);
		}
	    } while (pzomc < -1);
	  // F(EP) rejection
	  G4double XQC = 1.0+tau*(tau-2.0*cosTheta);
	  G4double AF = sqrt(XQC)*(1.0+tau*(tau-cosTheta)/XQC);
	  if (AF > 0) {
	    fpzmax = 1.0+AF*0.2;
	  }
	  else
	    {
	      fpzmax = 1.0-AF*0.2;
	    }
	  fpz = 1.0+AF*G4std::max(G4std::min(pzomc,0.2),-0.2);
	}while ((fpzmax*G4UniformRand())>fpz);
     
      //Energy of the scattered photon
      G4double T = pow(pzomc,2);
      G4double b1 = 1.0-T*tau*tau;
      G4double b2 = 1.0-T*tau*cosTheta;
      if (pzomc > 0.0)
	{
	  epsilon = (tau/b1)*(b2+sqrt(G4std::abs(b2*b2-b1*(1.0-T))));
	}
      else
	{
	  epsilon = (tau/b1)*(b2-sqrt(G4std::abs(b2*b2-b1*(1.0-T))));
	}
    }
  

  //G4cout << "Epsilon: " << epsilon << " cosTheta: " << cosTheta << G4endl;
  G4double sinTheta = sqrt(1-pow(cosTheta,2));
  G4double phi = twopi * G4UniformRand() ;
  G4double dirx = sinTheta * cos(phi);
  G4double diry = sinTheta * sin(phi);
  G4double dirz = cosTheta ;

  // Update G4VParticleChange for the scattered photon 
  
  G4ThreeVector photonDirection1(dirx,diry,dirz);
  photonDirection1.rotateUz(photonDirection0);
  aParticleChange.SetMomentumChange(photonDirection1) ;
  G4double photonEnergy1 = epsilon * photonEnergy0;   

  if (photonEnergy1 > 0.)
    {
      aParticleChange.SetEnergyChange(photonEnergy1) ;
    }
  else
    {    
      aParticleChange.SetEnergyChange(0.) ;
      aParticleChange.SetStatusChange(fStopAndKill);
    }

  //In teoria bisognerebbe generare anche la radiazione di fluorescenza

  // Kinematics of the scattered electron 

  G4double DiffEnergy = photonEnergy0*(1-epsilon);
  IonEnergy = (*((*IonizationEnergy)[Z-1]))[iosc];
  G4double eKineticEnergy = DiffEnergy - IonEnergy;
  G4double Q2 = pow(photonEnergy0,2)+photonEnergy1*(photonEnergy1-2.0*photonEnergy0*cosTheta);
  G4double cosThetaE; //scattering angle for the electron
  if (Q2 > 1.0e-12)
    {
      cosThetaE = (photonEnergy0-photonEnergy1*cosTheta)/sqrt(Q2);
    }
  else
    {
      cosThetaE = 1.0;
    }
  G4double sinThetaE = sqrt(1-pow(cosThetaE,2));
  // Generate the electron only if with large enough range w.r.t. cuts and safety

  G4double safety = aStep.GetPostStepPoint()->GetSafety();

  if (rangeTest->Escape(G4Electron::Electron(),material,eKineticEnergy,safety))
    {
      G4double xEl = sinThetaE * cos(phi+pi); 
      G4double yEl = sinThetaE * sin(phi+pi);
      G4double zEl = cosThetaE;
      G4ThreeVector eDirection(xEl,yEl,zEl); //electron direction
   
      G4DynamicParticle* electron = new G4DynamicParticle (G4Electron::Electron(),
							   eDirection,eKineticEnergy) ;
      aParticleChange.SetNumberOfSecondaries(1);
      aParticleChange.AddSecondary(electron);
      aParticleChange.SetLocalEnergyDeposit(0.); 
    }
  else
    {
      aParticleChange.SetNumberOfSecondaries(0);
      aParticleChange.SetLocalEnergyDeposit(eKineticEnergy);
    }

  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
}

G4bool G4PenelopeCompton::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() ); 
}

G4double G4PenelopeCompton::GetMeanFreePath(const G4Track& track, 
					     G4double previousStepSize, 
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
  G4String pathFile = pathString + "/penelope/comptondata.dat";
  G4std::ifstream file(pathFile);
  G4std::filebuf* lsdp = file.rdbuf();
  
  if (!(lsdp->is_open()))
    {
      G4String excep = "G4PenelopeCompton - data file " + pathFile + " not found!";
      G4Exception(excep);
    }
  file.close();
  
  //E' un algoritmo inguardabile...vedi se si puo' migliorare
  G4int k1,k2;
  G4double a1,a2;
  G4int Z=0;
  G4DataVector* f;
  G4DataVector* u;
  G4DataVector* j;
  for(Z=1;Z<93;Z++) {
    f = new G4DataVector;
    u = new G4DataVector;
    j = new G4DataVector;
    G4std::ifstream file(pathFile);
    k1=1;
    while (k1>0)
      {
	file >> k1 >> k2 >> a1 >> a2;
	if(k1 == Z)
	  {
	    f->push_back((G4double) k2);
	    u->push_back(a1);
	    j->push_back(a2);
	  }
      }
    file.close();
    IonizationEnergy->push_back(u);
    HartreeFunction->push_back(j);
    OccupationNumber->push_back(f);
  }
  //IonizationEnergy contiene Z=92 vettori di livelli di ionizzazione 
  //(*((*IonizationEnergy)[Z-1]))[i] contiene l'energia di ionizzazione dell'i=esimo
  //livello dell'elemento Z
};

G4double G4PenelopeCompton::CrossSection(G4double energy,G4int Z)
{
  G4double cs=0;
  if (energy<(5*MeV))
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
      G4double csl = 0.5*ki3*t0*t0+ki2*t0+ki1*log(t0)-(1.0/t0);
      G4int nosc = ((*OccupationNumber)[Z-1])->size();
      energyForIntegration=energy; 
      ZForIntegration = Z;
      for (G4int i=0;i<nosc;i++)
	{
	  G4double IonEnergy = (*((*IonizationEnergy)[Z-1]))[i];
	  G4double tau=(energy-IonEnergy)/energy;
	  if (tau < t0)
	    {
	      cs=pi*classic_electr_radius*classic_electr_radius*cs/(ki*ki3);
	      return cs;
		}
	  else
	    {
	      G4double csu = 0.5*ki3*tau*tau+ki2*tau+ki1*log(tau)-(1.0/tau);
	      G4int f = (G4int) (*((*OccupationNumber)[Z-1]))[i];
	      cs = cs + f*(csu-csl);
	    }
	}
      cs=pi*classic_electr_radius*classic_electr_radius*cs/(ki*ki3);
    }
  return cs;
}

  
G4double G4PenelopeCompton::DifferentialCrossSection(G4double cosTheta)
{
  const G4double k2 = sqrt(2.0);
  const G4double k1 = 1.0/k1;
  const G4double k12 = 0.5;
  G4double cdt1 = 1.0-cosTheta;
  G4double energy = energyForIntegration;
  G4int Z = ZForIntegration;
  G4double IonEnergy=0.0,Pzimax=0.0,XKN=0.0,DiffCS=0.0;
  G4double x=0.0,siap=0.0;
  G4double HarFunc=0.0;
  G4int OccupNb;
  //energy of Compton line;
  G4double EOEC = 1.0+(energy/electron_mass_c2)*cdt1; 
  G4double ECOE = 1.0/EOEC;
  //Incoherent scattering function (analytical profile)
  G4double sia = 0.0;
  G4int nosc = ((*OccupationNumber)[Z-1])->size();
  for (G4int i=0;i<nosc;i++){
    IonEnergy = (*((*IonizationEnergy)[Z-1]))[i];
    if (energy < IonEnergy) {
      XKN = EOEC+ECOE-1+cosTheta*cosTheta;
      DiffCS = pi*pow(classic_electr_radius,2)*pow(ECOE,2)*XKN*sia;
      return DiffCS;
    }
    else
      {
	G4double aux = energy * (energy-IonEnergy)*cdt1;
	Pzimax = (aux - electron_mass_c2*IonEnergy)/(electron_mass_c2*sqrt(2*aux+pow(IonEnergy,2)));
	HarFunc = (*((*HartreeFunction)[Z-1]))[i];
	OccupNb = (G4int) (*((*OccupationNumber)[Z-1]))[i];
	x = HarFunc*Pzimax;
	if (x > 0) 
	  {
	    siap = 1.0-0.5*exp(k12-pow((k1+k2*x),2));
	  }
	else
	  {
	    siap = 0.5*exp(k12-pow((k1-k2*x),2));
	  }
	sia = sia + OccupNb*siap; //sum of all contributions;
      }
  }
  XKN = EOEC+ECOE-1+cosTheta*cosTheta;
  DiffCS = pi*pow(classic_electr_radius,2)*pow(ECOE,2)*XKN*sia;
  return DiffCS;
}

G4int G4PenelopeCompton::SelectRandomAtomForCompton(G4Material* material,G4double energy) const
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
  for (G4int i=0;i<nElements;i++)
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
  
