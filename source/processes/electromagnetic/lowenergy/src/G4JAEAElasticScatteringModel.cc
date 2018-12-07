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
M. Omer and R. Hajima  on   17 October 2016
contact:
omer.mohamed@jaea.go.jp and hajima.ryoichi@qst.go.jp
Publication Information:
1- M. Omer, R. Hajima, Including Delbrück scattering in Geant4,
Nucl. Instrum. Methods Phys. Res. Sect. B, vol. 405, 2017, pp. 43-49.,
https://doi.org/10.1016/j.nimb.2017.05.028
2- M. Omer, R. Hajima, Geant4 physics process for elastic scattering of gamma-rays,
JAEA Technical Report 2018-007, 2018.
https://doi.org/10.11484/jaea-data-code-2018-007
*/
//         on base of G4LivermoreRayleighModel


#include "G4JAEAElasticScatteringModel.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4JAEAElasticScatteringModel::maxZ = 99;
G4LPhysicsFreeVector* G4JAEAElasticScatteringModel::dataCS[] ;
//Initialising an array to hold all elastic scattering data.
G4double Diff_CS_data[100][183][300];

G4JAEAElasticScatteringModel::G4JAEAElasticScatteringModel()
  :G4VEmModel("G4JAEAElasticScatteringModel"),isInitialised(false)
{
  fParticleChange = 0;
//Low energy limit for G4JAEAElasticScatteringModel process.
  lowEnergyLimit  = 10 * keV;

  verboseLevel= 0;
  // Verbosity scale for debugging purposes:
  // 0 = nothing
  // 1 = calculation of cross sections, file openings...
  // 2 = entering in methods

  if(verboseLevel > 0)
  {
    G4cout << "G4JAEAElasticScatteringModel is constructed " << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4JAEAElasticScatteringModel::~G4JAEAElasticScatteringModel()
{
  if(IsMaster()) {
    for(G4int i=0; i<maxZ; ++i) {
      if(dataCS[i]) {
	delete dataCS[i];
	dataCS[i] = 0;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4JAEAElasticScatteringModel::Initialise(const G4ParticleDefinition* particle,
					  const G4DataVector& cuts)
{
  if (verboseLevel > 1)
  {
    G4cout << "Calling Initialise() of G4JAEAElasticScatteringModel." << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / eV << " eV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }

  if(IsMaster()) {

    // Initialise element selector
    InitialiseElementSelectors(particle, cuts);

    // Access to elements
    char* path = getenv("G4LEDATA");
    G4ProductionCutsTable* theCoupleTable =
      G4ProductionCutsTable::GetProductionCutsTable();
    G4int numOfCouples = theCoupleTable->GetTableSize();

    for(G4int i=0; i<numOfCouples; ++i)
      {
	const G4MaterialCutsCouple* couple =
	  theCoupleTable->GetMaterialCutsCouple(i);
	const G4Material* material = couple->GetMaterial();
	const G4ElementVector* theElementVector = material->GetElementVector();
	G4int nelm = material->GetNumberOfElements();

	for (G4int j=0; j<nelm; ++j)
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

void G4JAEAElasticScatteringModel::InitialiseLocal(const G4ParticleDefinition*,
					       G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4JAEAElasticScatteringModel::ReadData(size_t Z, const char* path)
{

 if (verboseLevel > 1)
  {
    G4cout << "Calling ReadData() of G4JAEAElasticScatteringModel"
	   << G4endl;
  }

  if(dataCS[Z]) { return; }

  const char* datadir = path;

  if(!datadir)
  {
    datadir = getenv("G4LEDATA");
    if(!datadir)
    {
      G4Exception("G4JAEAElasticScatteringModel::ReadData()","em0006",
		  FatalException,
		  "Environment variable G4LEDATA not defined");
      return;
    }
  }

/* Reading all data in the form of 183 * 300 array.
The first row is the energy, and the second row is the total cross section.
Rows from the 3rd to the 183rd are the differential cross section with an angular resolution of 1 degree.
*/
G4double ESdata[183][300];

std::ostringstream ostCS;
ostCS << datadir << "/JAEAESData/cs_Z_" << Z <<".dat";
std::ifstream alldata(ostCS.str().c_str());
if(!alldata.is_open())
{
  G4ExceptionDescription ed;
  ed << "G4JAEAElasticScattering Model data file <" << ostCS.str().c_str()
     << "> is not opened!" << G4endl;
  G4Exception("Elastic Scattering::ReadData()","em0003",FatalException,
  ed,"G4LEDATA version should be G4EMLOW6.27 or later. Elastic Scattering Data are not loaded");
  return;
}
else
{
  if(verboseLevel > 3) {
    G4cout << "File " << ostCS.str()
	   << " is opened by G4JAEAElasticScatteringModel" << G4endl;
  }
  }
while (!alldata.eof())
{
	for (int i=0; i<183;i++)
	{
		for (int j=0; j<300; j++)
		{
		alldata >> ESdata[i][j];
		Diff_CS_data[Z][i][j]=ESdata[i][j];
		}
	}
	if (!alldata) break;
}

/*
Writing the total cross section data to a G4LPhysicsFreeVector.
This provides an interpolation of the Energy-Total Cross Section data.
*/

      dataCS[Z] = new G4LPhysicsFreeVector(300,0.01,3.);
//Note that the total cross section and energy are converted to the internal units.
      for (int i=0;i<300;i++)
      dataCS[Z]->PutValue(i,Diff_CS_data[Z][0][i]*1e-3,Diff_CS_data[Z][1][i]*1e-22);

  // Activation of spline interpolation
    dataCS[Z] ->SetSpline(true);



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4JAEAElasticScatteringModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{

  if (verboseLevel > 1)
  {
    G4cout << "G4JAEAElasticScatteringModel::ComputeCrossSectionPerAtom()"
	   << G4endl;
  }

  if(GammaEnergy < lowEnergyLimit) { return 0.0; }

  G4double xs = 0.0;

  G4int intZ = G4lrint(Z);

  if(intZ < 1 || intZ > maxZ) { return xs; }

  G4LPhysicsFreeVector* pv = dataCS[intZ];

  // if element was not initialised
  // do initialisation safely for MT mode
  if(!pv) {
    InitialiseForElement(0, intZ);
    pv = dataCS[intZ];
    if(!pv) { return xs; }
  }

  G4int n = pv->GetVectorLength() - 1;

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

void G4JAEAElasticScatteringModel::SampleSecondaries(
                          std::vector<G4DynamicParticle*>*,
			  const G4MaterialCutsCouple* couple,
			  const G4DynamicParticle* aDynamicGamma,
			  G4double, G4double)
{
  if (verboseLevel > 1) {
    G4cout << "Calling SampleSecondaries() of G4JAEAElasticScatteringModel"
	   << G4endl;
  }
  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();

  // Absorption of low-energy gamma
  if (photonEnergy0 <= lowEnergyLimit)
    {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0);
      return ;
    }


  // Select randomly one element in the current material
  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple,particle,photonEnergy0);
  G4int Z = G4lrint(elm->GetZ());


//Select the angular distribution depending on the photon energy
	G4double *whichdistribution = lower_bound(Diff_CS_data[Z][0],Diff_CS_data[Z][0]+300,photonEnergy0*1000.);
	int index = max(0,(int)(whichdistribution-Diff_CS_data[Z][0]-1));

//Rounding up to half the energy-grid separation (5 keV)
	if (photonEnergy0*1000>=0.5*(Diff_CS_data[Z][0][index]+Diff_CS_data[Z][0][index+1]))
			index++;
/*
Getting the normalized probablity distrbution function and
normalization factor to create the probability distribution function
*/
	G4double normdist=0;
	for (int i=0;i<=180;i++)
		{
      		distribution[i]=Diff_CS_data[Z][i+2][index];
		normdist = normdist + distribution[i];
		}
//Create the cummulative distribution function (cdf)
	for (int i =0;i<=180;i++) pdf[i]=distribution[i]/normdist;
	cdf[0]=0;
	G4double cdfsum =0;
	for (int i=0; i<=180;i++)
		{
		cdfsum=cdfsum+pdf[i];
		cdf[i]=cdfsum;
		}
//Sampling the polar angle by inverse transform uing cdf.
    G4double r = G4UniformRand();
    G4double *cdfptr=lower_bound(cdf,cdf+181,r);
    int cdfindex = (int)(cdfptr-cdf-1);
    G4double cdfinv = (r-cdf[cdfindex])/(cdf[cdfindex+1]-cdf[cdfindex]);
    G4double theta = (cdfindex+cdfinv)/180.;
//polar is now ready
    theta = theta*CLHEP::pi;


/* Alternative sampling using CLHEP functions
    CLHEP::RandGeneral GenDistTheta(distribution,181);
    G4double theta = CLHEP::pi*GenDistTheta.shoot();
   theta =theta*CLHEP::pi;		//polar is now ready
*/

//Azimuth is uniformally distributed
G4double phi  = CLHEP::twopi*G4UniformRand();

G4ThreeVector finaldirection(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
finaldirection.rotateUz(aDynamicGamma->GetMomentumDirection());
//Sampling the Final State
		fParticleChange->ProposeMomentumDirection(finaldirection);
		fParticleChange->SetProposedKineticEnergy(photonEnergy0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AutoLock.hh"
namespace { G4Mutex G4JAEAElasticScatteringModelMutex = G4MUTEX_INITIALIZER; }

void
G4JAEAElasticScatteringModel::InitialiseForElement(const G4ParticleDefinition*,
					       G4int Z)
{
  G4AutoLock l(&G4JAEAElasticScatteringModelMutex);
  if(!dataCS[Z]) { ReadData(Z); }
  l.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
