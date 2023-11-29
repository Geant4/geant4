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
/*
  Authors:
  M. Omer and R. Hajima  on 15 November 2019
  contact:
  omer.mohamed@jaea.go.jp and hajima.ryoichi@qst.go.jp
  Publication Information:
  1- M. Omer, R. Hajima, Validating polarization effects in gamma-rays elastic scattering by Monte
  Carlo simulation, New J. Phys., vol. 21, 2019, pp. 113006 (1-10),
  https://doi.org/10.1088/1367-2630/ab4d8a
*/

#include "G4JAEAPolarizedElasticScatteringModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"
namespace { G4Mutex G4JAEAPolarizedElasticScatteringModelMutex = G4MUTEX_INITIALIZER; }
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsFreeVector* G4JAEAPolarizedElasticScatteringModel::dataCS[] = {nullptr};
G4DataVector* G4JAEAPolarizedElasticScatteringModel::Polarized_ES_Data[] = {nullptr};

G4JAEAPolarizedElasticScatteringModel::G4JAEAPolarizedElasticScatteringModel()
  :G4VEmModel("G4JAEAPolarizedElasticScatteringModel"),isInitialised(false)
{
  fParticleChange = 0;
  lowEnergyLimit  = 100 * keV;	 //low energy limit for JAEAElasticScattering cross section data
  fLinearPolarizationSensitvity1=1;
  fLinearPolarizationSensitvity2=1;
  fCircularPolarizationSensitvity=1;

  verboseLevel= 0;
  // Verbosity scale for debugging purposes:
  // 0 = nothing
  // 1 = calculation of cross sections, file openings...
  // 2 = entering in methods

  if(verboseLevel > 0)
    {
      G4cout << "G4JAEAPolarizedElasticScatteringModel is constructed " << G4endl;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4JAEAPolarizedElasticScatteringModel::~G4JAEAPolarizedElasticScatteringModel()
{
  if(IsMaster()) {
    for(G4int i=0; i<=maxZ; ++i) {
      if(dataCS[i]) {
	delete dataCS[i];
	dataCS[i] = nullptr;
      }
      if (Polarized_ES_Data[i]){
	delete Polarized_ES_Data[i];
	Polarized_ES_Data[i] = nullptr;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4JAEAPolarizedElasticScatteringModel::Initialise(const G4ParticleDefinition* particle,
						       const G4DataVector& cuts)
{
  if (verboseLevel > 1)
    {
      G4cout << "Calling Initialise() of G4JAEAPolarizedElasticScatteringModel." << G4endl
	     << "Energy range: "
	     << LowEnergyLimit() / eV << " eV - "
	     << HighEnergyLimit() / GeV << " GeV"
	     << G4endl;
    }

  if(IsMaster()) {

    // Initialise element selector
    InitialiseElementSelectors(particle, cuts);

    // Access to elements
    const char* path = G4FindDataDir("G4LEDATA");
    G4ProductionCutsTable* theCoupleTable =
      G4ProductionCutsTable::GetProductionCutsTable();
    G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

    for(G4int i=0; i<numOfCouples; ++i)
      {
	const G4MaterialCutsCouple* couple =
	  theCoupleTable->GetMaterialCutsCouple(i);
	const G4Material* material = couple->GetMaterial();
	const G4ElementVector* theElementVector = material->GetElementVector();
	std::size_t nelm = material->GetNumberOfElements();

	for (std::size_t j=0; j<nelm; ++j)
	  {
	    G4int Z = G4lrint((*theElementVector)[j]->GetZ());
	    if(Z < 1)          { Z = 1; }
	    else if(Z > maxZ)  { Z = maxZ; }
	    if( (!dataCS[Z]) ) { ReadData(Z, path); }
	  }
      }
  }

  if(isInitialised) { return; }
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4JAEAPolarizedElasticScatteringModel::InitialiseLocal(const G4ParticleDefinition*,
							    G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4JAEAPolarizedElasticScatteringModel::ReadData(std::size_t Z, const char* path)
{
  if (verboseLevel > 1)
    {
      G4cout << "Calling ReadData() of G4JAEAPolarizedElasticScatteringModel"
	     << G4endl;
    }

  if(dataCS[Z]) { return; }

  const char* datadir = path;
  if(!datadir)
    {
      datadir = G4FindDataDir("G4LEDATA");
      if(!datadir)
	{
	  G4Exception("G4JAEAPolarizedElasticScatteringModel::ReadData()","em0006",
		      FatalException,
		      "Environment variable G4LEDATA not defined");
	  return;
	}
    }
  
  std::ostringstream ostCS;
  ostCS << datadir << "/JAEAESData/amp_Z_" << Z ;
  std::ifstream ES_Data_Buffer(ostCS.str().c_str(),ios::binary);
  if( !ES_Data_Buffer.is_open() )
    {
      G4ExceptionDescription ed;
      ed << "G4JAEAPolarizedElasticScattering Model data file <" << ostCS.str().c_str()
	 << "> is not opened!" << G4endl;
      G4Exception("G4JAEAPolarizedElasticScatteringModel::ReadData()","em0003",FatalException,
		  ed,"G4LEDATA version should be G4EMLOW7.11 or later. Polarized Elastic Scattering Data are not loaded");
      return;
    }
  else
    {
      if(verboseLevel > 3) {
	G4cout << "File " << ostCS.str()
	       << " is opened by G4JAEAPolarizedElasticScatteringModel" << G4endl;
      }
    }
  
  
  if (!Polarized_ES_Data[Z])
    Polarized_ES_Data[Z] = new G4DataVector();
  
  G4float buffer_var;
  while (ES_Data_Buffer.read(reinterpret_cast<char*>(&buffer_var),sizeof(float)))
    {
      Polarized_ES_Data[Z]->push_back(buffer_var);
    }
  
  dataCS[Z] = new G4PhysicsFreeVector(300,0.01,3.,/*spline=*/true);
  
  for (G4int i=0;i<300;++i)
    dataCS[Z]->PutValues(i,10.*i*1e-3,Polarized_ES_Data[Z]->at(i)*1e-22);
 
  // Activation of spline interpolation
  dataCS[Z] ->FillSecondDerivatives();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4JAEAPolarizedElasticScatteringModel::ComputeCrossSectionPerAtom(
									   const G4ParticleDefinition*,
									   G4double GammaEnergy,
									   G4double Z, G4double,
									   G4double, G4double)
{
  //Select the energy-grid point closest to the photon energy
  //		G4double *whichenergy = lower_bound(ESdata[0],ESdata[0]+300,GammaEnergy);
  //		int energyindex = max(0,(int)(whichenergy-ESdata[0]-1));
  
  if (verboseLevel > 1)
    {
      G4cout << "G4JAEAPolarizedElasticScatteringModel::ComputeCrossSectionPerAtom()"
	     << G4endl;
    }
  
  if(GammaEnergy < lowEnergyLimit) { return 0.0; }
  
  G4double xs = 0.0;
  
  G4int intZ = G4lrint(Z);
  
  if(intZ < 1 || intZ > maxZ) { return xs; }

  G4PhysicsFreeVector* pv = dataCS[intZ];

  // if element was not initialised
  // do initialisation safely for MT mode
  if(!pv) {
    InitialiseForElement(0, intZ);
    pv = dataCS[intZ];
    if(!pv) { return xs; }
  }

  std::size_t n = pv->GetVectorLength() - 1;

  G4double e = GammaEnergy;
  if(e >= pv->Energy(n)) {
    xs = (*pv)[n];
  } else if(e >= pv->Energy(0)) {
    xs = pv->Value(e);
  }

  if(verboseLevel > 0)
    {
      G4cout  <<  "****** DEBUG: tcs value for Z=" << Z << " at energy (MeV)="
	      << e << G4endl;
      G4cout  <<  "  cs (Geant4 internal unit)=" << xs << G4endl;
      G4cout  <<  "    -> first E*E*cs value in CS data file (iu) =" << (*pv)[0]
	      << G4endl;
      G4cout  <<  "    -> last  E*E*cs value in CS data file (iu) =" << (*pv)[n]
	      << G4endl;
      G4cout  <<  "*********************************************************"
	      << G4endl;
    }

  return (xs);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4JAEAPolarizedElasticScatteringModel::SampleSecondaries(
							      std::vector<G4DynamicParticle*>*,
							      const G4MaterialCutsCouple* couple,
							      const G4DynamicParticle* aDynamicGamma,
							      G4double, G4double)
{
  if (verboseLevel > 1) {

    G4cout << "Calling SampleSecondaries() of G4JAEAPolarizedElasticScatteringModel."
	   << G4endl;
  }
  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();

  // absorption of low-energy gamma
  if (photonEnergy0 <= lowEnergyLimit)
    {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0);
      return ;
    }

  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple,particle,photonEnergy0);
  G4int Z = G4lrint(elm->GetZ());

  //Getting the corresponding distrbution
  G4int energyindex=round(100*photonEnergy0)-1;
  G4double a1=0, a2=0, a3=0,a4=0;
  for (G4int i=0;i<=180;++i)
    {
      a1=Polarized_ES_Data[Z]->at(4*i+300+181*4*(energyindex));
      a2=Polarized_ES_Data[Z]->at(4*i+1+300+181*4*(energyindex));
      a3=Polarized_ES_Data[Z]->at(4*i+2+300+181*4*(energyindex));
      a4=Polarized_ES_Data[Z]->at(4*i+3+300+181*4*(energyindex));
      distribution[i]=a1*a1+a2*a2+a3*a3+a4*a4;
    }

  CLHEP::RandGeneral GenThetaDist(distribution,180);
  //Intial sampling of the scattering angle. To be updated for the circular polarization
  G4double theta = CLHEP::pi*GenThetaDist.shoot();
  //G4double theta =45.*CLHEP::pi/180.;
  //Theta is in degree to call scattering amplitudes
  G4int theta_in_degree =round(theta*180./CLHEP::pi);

  //theta_in_degree=45;

  G4double am1=0,am2=0,am3=0,am4=0,aparaSquare=0,aperpSquare=0,
    apara_aper_Asterisk=0,img_apara_aper_Asterisk=0;
  am1=Polarized_ES_Data[Z]->at(4*theta_in_degree+300+181*4*(energyindex));
  am2=Polarized_ES_Data[Z]->at(4*theta_in_degree+1+300+181*4*(energyindex));
  am3=Polarized_ES_Data[Z]->at(4*theta_in_degree+2+300+181*4*(energyindex));
  am4=Polarized_ES_Data[Z]->at(4*theta_in_degree+3+300+181*4*(energyindex));
  aparaSquare=am1*am1+am2*am2;
  aperpSquare=am3*am3+am4*am4;
  apara_aper_Asterisk=2*a1*a3+2*a2*a4;
  img_apara_aper_Asterisk=2*a1*a4-2*a2*a3;

  G4ThreeVector Direction_Unpolarized(0.,0.,0.);
  G4ThreeVector Direction_Linear1(0.,0.,0.);
  G4ThreeVector Direction_Linear2(0.,0.,0.);
  G4ThreeVector Direction_Circular(0.,0.,0.);
  G4ThreeVector Polarization_Unpolarized(0.,0.,0.);
  G4ThreeVector Polarization_Linear1(0.,0.,0.);
  G4ThreeVector Polarization_Linear2(0.,0.,0.);
  G4ThreeVector Polarization_Circular(0.,0.,0.);

  //Stokes parameters for the incoming and outgoing photon
  G4double Xi1=0, Xi2=0, Xi3=0, Xi1_Prime=0,Xi2_Prime=0,Xi3_Prime=0;
    
  //Getting the Stokes parameters for the incoming photon
  G4ThreeVector gammaPolarization0 = aDynamicGamma->GetPolarization();
  Xi1=gammaPolarization0.x();
  Xi2=gammaPolarization0.y();
  Xi3=gammaPolarization0.z();
    
  //Polarization vector must be unit vector (5% tolerance)    
  if ((gammaPolarization0.mag())>1.05 || (Xi1*Xi1>1.05) || (Xi2*Xi2>1.05) || (Xi3*Xi3>1.05))
    {
      G4Exception("G4JAEAPolarizedElasticScatteringModel::SampleSecondaries()","em1006",
		  JustWarning,
		  "WARNING: G4JAEAPolarizedElasticScatteringModel is only compatible with a unit polarization vector.");
      return;
    }
  //Unpolarized gamma rays
  if (Xi1==0 && Xi2==0 && Xi3==0)
    {
      G4double Phi_Unpolarized=0;
      if (fLinearPolarizationSensitvity1)
	Phi_Unpolarized=GeneratePolarizedPhi(aparaSquare,aperpSquare,0.);
      else
	Phi_Unpolarized=CLHEP::twopi*G4UniformRand();
      Direction_Unpolarized.setX(sin(theta)*cos(Phi_Unpolarized));
      Direction_Unpolarized.setY(sin(theta)*sin(Phi_Unpolarized));
      Direction_Unpolarized.setZ(cos(theta));
      Direction_Unpolarized.rotateUz(aDynamicGamma->GetMomentumDirection());
      Xi1_Prime=(aparaSquare-aperpSquare)/(aparaSquare+aperpSquare);
      Polarization_Unpolarized.setX(Xi1_Prime);
      Polarization_Unpolarized.setY(0.);
      Polarization_Unpolarized.setZ(0.);
      fParticleChange->ProposeMomentumDirection(Direction_Unpolarized);
      fParticleChange->ProposePolarization(Polarization_Unpolarized);
      return;
    }

  //Linear polarization defined by first Stokes parameter
  G4double InitialAzimuth=aDynamicGamma->GetMomentumDirection().phi();
  if(InitialAzimuth<0) InitialAzimuth=InitialAzimuth+CLHEP::twopi;
    
  G4double Phi_Linear1=0.;

  Phi_Linear1 = GeneratePolarizedPhi(aparaSquare+aperpSquare+Xi1*(aparaSquare-aperpSquare),
				     aparaSquare+aperpSquare-Xi1*(aparaSquare-aperpSquare),InitialAzimuth);

  Xi1_Prime=((aparaSquare-aperpSquare)+Xi1*(aparaSquare+aperpSquare)*cos(2*Phi_Linear1))/
    ((aparaSquare+aperpSquare)+Xi1*(aparaSquare-aperpSquare)*cos(2*Phi_Linear1));
  Xi2_Prime=(-Xi1*apara_aper_Asterisk*sin(2*Phi_Linear1))/
    ((aparaSquare+aperpSquare)+Xi1*(aparaSquare-aperpSquare)*cos(2*Phi_Linear1));
  Xi3_Prime=(-Xi1*img_apara_aper_Asterisk*sin(2*Phi_Linear1))/
    ((aparaSquare+aperpSquare)+Xi1*(aparaSquare-aperpSquare)*cos(2*Phi_Linear1));
  //Store momentum direction and po;arization
  Direction_Linear1.setX(sin(theta)*cos(Phi_Linear1));
  Direction_Linear1.setY(sin(theta)*sin(Phi_Linear1));
  Direction_Linear1.setZ(cos(theta));
  Polarization_Linear1.setX(Xi1_Prime);
  Polarization_Linear1.setY(Xi2_Prime);
  Polarization_Linear1.setZ(Xi3_Prime);
    
  //Set scattered photon polarization sensitivity
  Xi1_Prime=Xi1_Prime*fLinearPolarizationSensitvity1;
  Xi2_Prime=Xi2_Prime*fLinearPolarizationSensitvity2;
  Xi3_Prime=Xi3_Prime*fCircularPolarizationSensitvity;
    
  G4double dsigmaL1=0.0;
  if(abs(Xi1)>0.0) dsigmaL1=0.25*((aparaSquare+aperpSquare)*
				  (1+Xi1*Xi1_Prime*cos(2*Phi_Linear1))+
				  (aparaSquare-aperpSquare)*(Xi1*cos(2*Phi_Linear1)+Xi1_Prime)
				  -Xi1*Xi2_Prime*apara_aper_Asterisk*sin(2*Phi_Linear1)-
				  Xi1*Xi3_Prime*img_apara_aper_Asterisk*sin(2*Phi_Linear1));

  //Linear polarization defined by second Stokes parameter
  //G4double IntialAzimuth=aDynamicGamma->GetMomentumDirection().phi();
  G4double Phi_Linear2=0.;

  InitialAzimuth=InitialAzimuth-CLHEP::pi/4.;
  if(InitialAzimuth<0) InitialAzimuth=InitialAzimuth+CLHEP::twopi;

  Phi_Linear2 = GeneratePolarizedPhi(aparaSquare+aperpSquare+Xi1*(aparaSquare-aperpSquare)
				     ,aparaSquare+aperpSquare-Xi1*(aparaSquare-aperpSquare),InitialAzimuth);
    
  Xi1_Prime=((aparaSquare-aperpSquare)+Xi2*(aparaSquare+aperpSquare)*sin(2*Phi_Linear2))/
    ((aparaSquare+aperpSquare)+Xi2*(aparaSquare-aperpSquare)*sin(2*Phi_Linear2));
  Xi2_Prime=(Xi2*apara_aper_Asterisk*cos(2*Phi_Linear2))/
    ((aparaSquare+aperpSquare)+Xi2*(aparaSquare-aperpSquare)*sin(2*Phi_Linear2));
  Xi3_Prime=(Xi2*img_apara_aper_Asterisk*cos(2*Phi_Linear2))/
    ((aparaSquare+aperpSquare)+Xi2*(aparaSquare-aperpSquare)*sin(2*Phi_Linear2));
  //Store momentum direction and polarization
  Direction_Linear2.setX(sin(theta)*cos(Phi_Linear2));
  Direction_Linear2.setY(sin(theta)*sin(Phi_Linear2));
  Direction_Linear2.setZ(cos(theta));
  Polarization_Linear2.setX(Xi1_Prime);
  Polarization_Linear2.setY(Xi2_Prime);
  Polarization_Linear2.setZ(Xi3_Prime);

  //Set scattered photon polarization sensitivity
  Xi1_Prime=Xi1_Prime*fLinearPolarizationSensitvity1;
  Xi2_Prime=Xi2_Prime*fLinearPolarizationSensitvity2;
  Xi3_Prime=Xi3_Prime*fCircularPolarizationSensitvity;

  G4double dsigmaL2=0.0;
  if(abs(Xi2)>0.0)
    dsigmaL2=0.25*((aparaSquare+aperpSquare)*(1+Xi2*Xi1_Prime*sin(2*Phi_Linear2))+
		   (aparaSquare-aperpSquare)*(Xi2*sin(2*Phi_Linear2)+Xi1_Prime)
		   +Xi2*Xi2_Prime*apara_aper_Asterisk*cos(2*Phi_Linear2)-
		   Xi2*Xi3_Prime*img_apara_aper_Asterisk*cos(2*Phi_Linear2));

  //Circular polarization
  G4double Phi_Circular = CLHEP::twopi*G4UniformRand();
  G4double Theta_Circular = 0.;
    
  Xi1_Prime=(aparaSquare-aperpSquare)/(aparaSquare+aperpSquare);
  Xi2_Prime=(-Xi3*img_apara_aper_Asterisk)/(aparaSquare+aperpSquare);
  Xi3_Prime=(Xi3*apara_aper_Asterisk)/(aparaSquare+aperpSquare);
    
  Polarization_Circular.setX(Xi1_Prime);
  Polarization_Circular.setY(Xi2_Prime);
  Polarization_Circular.setZ(Xi3_Prime);
    
  //Set scattered photon polarization sensitivity
  Xi1_Prime=Xi1_Prime*fLinearPolarizationSensitvity1;
  Xi2_Prime=Xi2_Prime*fLinearPolarizationSensitvity2;
  Xi3_Prime=Xi3_Prime*fCircularPolarizationSensitvity;
    
  G4double dsigmaC=0.0;
  if(abs(Xi3)>0.0)
    dsigmaC=0.25*(aparaSquare+aperpSquare+Xi1_Prime*(aparaSquare-aperpSquare)-
		  Xi3*Xi2_Prime*img_apara_aper_Asterisk
		  +Xi3*Xi3_Prime*apara_aper_Asterisk);

  if (abs(Xi3)==0.0 && abs(Xi1_Prime)==0.0)
    {
      Direction_Circular.setX(sin(theta)*cos(Phi_Circular));
      Direction_Circular.setY(sin(theta)*sin(Phi_Circular));
      Direction_Circular.setZ(cos(theta));
    }
  else
    {
      G4double c1=0, c2=0, c3=0,c4=0;
      for (G4int i=0;i<=180;++i)
	{
	  c1=Polarized_ES_Data[Z]->at(4*i+300+181*4*(energyindex));
	  c2=Polarized_ES_Data[Z]->at(4*i+1+300+181*4*(energyindex));
	  c3=Polarized_ES_Data[Z]->at(4*i+2+300+181*4*(energyindex));
	  c4=Polarized_ES_Data[Z]->at(4*i+3+300+181*4*(energyindex));
	  cdistribution[i]=0.25*((c1*c1+c2*c2+c3*c3+c4*c4)+
				 Xi1_Prime*(c1*c1+c2*c2-c3*c3-c4*c4)-
				 Xi3*Xi2_Prime*(2*c1*c4-2*c2*c3)
				 +Xi3*Xi3_Prime*(2*c1*c4-2*c2*c3));
	}
      CLHEP::RandGeneral GenTheta_Circ_Dist(cdistribution,180);
      Theta_Circular=CLHEP::pi*GenTheta_Circ_Dist.shoot();
      Direction_Circular.setX(sin(Theta_Circular)*cos(Phi_Circular));
      Direction_Circular.setY(sin(Theta_Circular)*sin(Phi_Circular));
      Direction_Circular.setZ(cos(Theta_Circular));
    }

  // Sampling scattered photon direction based on asymmetry arising from polarization mixing
  G4double totalSigma= dsigmaL1+dsigmaL2+dsigmaC;
  G4double prob1=dsigmaL1/totalSigma;
  G4double prob2=dsigmaL2/totalSigma;
  G4double probc=1-(prob1+prob2);
    
  //Check the Probability of polarization mixing
  if (abs(probc - dsigmaC/totalSigma)>=0.0001)
    {
      G4Exception("G4JAEAPolarizedElasticScatteringModel::SampleSecondaries()","em1007",
		  JustWarning,
		  "WARNING: Polarization mixing might be incorrect.");
    }
    
  // Generate outgoing photon direction
  G4ThreeVector finaldirection(0.0,0.0,0.0);
  G4ThreeVector outcomingPhotonPolarization(0.0,0.0,0.0);
    
  //Polarization mixing
  G4double polmix=G4UniformRand();
  if (polmix<=prob1)
    {
      finaldirection.setX(Direction_Linear1.x());
      finaldirection.setY(Direction_Linear1.y());
      finaldirection.setZ(Direction_Linear1.z());
      outcomingPhotonPolarization.setX(Polarization_Linear1.x());
      outcomingPhotonPolarization.setY(Polarization_Linear1.y());
      outcomingPhotonPolarization.setZ(Polarization_Linear1.z());
    }
  else if ((polmix>prob1) && (polmix<=prob1+prob2))
    {
      finaldirection.setX(Direction_Linear2.x());
      finaldirection.setY(Direction_Linear2.y());
      finaldirection.setZ(Direction_Linear2.z());
      outcomingPhotonPolarization.setX(Polarization_Linear2.x());
      outcomingPhotonPolarization.setY(Polarization_Linear2.y());
      outcomingPhotonPolarization.setZ(Polarization_Linear2.z());
    }
  else if (polmix>prob1+prob2)
    {
      finaldirection.setX(Direction_Circular.x());
      finaldirection.setY(Direction_Circular.y());
      finaldirection.setZ(Direction_Circular.z());
      outcomingPhotonPolarization.setX(Polarization_Circular.x());
      outcomingPhotonPolarization.setY(Polarization_Circular.y());
      outcomingPhotonPolarization.setZ(Polarization_Circular.z());
    }

  //Sampling the Final State
  finaldirection.rotateUz(aDynamicGamma->GetMomentumDirection());
  fParticleChange->ProposeMomentumDirection(finaldirection);
  fParticleChange->SetProposedKineticEnergy(photonEnergy0);
  fParticleChange->ProposePolarization(outcomingPhotonPolarization);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4JAEAPolarizedElasticScatteringModel::GeneratePolarizedPhi(G4double Sigma_para,
								     G4double Sigma_perp, 
								     G4double initial_Pol_Plane)
{
  G4double phi;
  G4double phiProbability;
  G4double Probability=Sigma_perp/(Sigma_para+Sigma_perp);
  if (Probability<=G4UniformRand())
    {
      do
	{
	  phi = CLHEP::twopi * G4UniformRand();
	  phiProbability = cos(phi+initial_Pol_Plane)*cos(phi+initial_Pol_Plane);
	}
      while (phiProbability < G4UniformRand());

    }
  else
    {
      do
	{
	  phi = CLHEP::twopi * G4UniformRand();
	  phiProbability = sin(phi+initial_Pol_Plane)*sin(phi+initial_Pol_Plane);
	}
      while (phiProbability < G4UniformRand());
    }
  return phi;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
G4JAEAPolarizedElasticScatteringModel::InitialiseForElement(const G4ParticleDefinition*,
							    G4int Z)
{
  G4AutoLock l(&G4JAEAPolarizedElasticScatteringModelMutex);
  //  G4cout << "G4JAEAPolarizedElasticScatteringModel::InitialiseForElement Z= "
  //   << Z << G4endl;
  if(!dataCS[Z]) { ReadData(Z); }
  l.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
