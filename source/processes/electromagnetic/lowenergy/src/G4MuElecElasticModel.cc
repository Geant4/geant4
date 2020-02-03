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
// G4MuElecElasticModel.cc, 2011/08/29 A.Valentin, M. Raine
//
// Based on the following publications
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for electrons in Si,
//	    NIM B, vol. 288, pp. 66 - 73, 2012.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "G4MuElecElasticModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuElecElasticModel::G4MuElecElasticModel(const G4ParticleDefinition*,
                                             const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The name of the class G4MuElecElasticModel is changed to G4MicroElecElasticModel. " << G4endl;
   G4cout << "   The obsolete class will be REMOVED with the next release of Geant4. " << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;

  nistSi = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

  killBelowEnergy = 16.7 * eV; // Minimum e- energy for energy loss by excitation
  lowEnergyLimit = 0 * eV;
  lowEnergyLimitOfModel = 5 * eV; // The model lower energy is 5 eV
  highEnergyLimit = 100. * MeV;
  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);

  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if( verboseLevel>0 )
  {
    G4cout << "MuElec Elastic model is constructed " << G4endl
           << "Energy range: "
           << lowEnergyLimit / eV << " eV - "
           << highEnergyLimit / keV << " keV"
           << G4endl;
  }
  fParticleChangeForGamma = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuElecElasticModel::~G4MuElecElasticModel()
{
  // For total cross section

  std::map< G4String,G4MuElecCrossSectionDataSet*,std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4MuElecCrossSectionDataSet* table = pos->second;
    delete table;
  }

   // For final state

   eVecm.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuElecElasticModel::Initialise(const G4ParticleDefinition* /*particle*/,
                                       const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
    G4cout << "Calling G4MuElecElasticModel::Initialise()" << G4endl;

  // Energy limits

  if (LowEnergyLimit() < lowEnergyLimit)
  {
    G4cout << "G4MuElecElasticModel: low energy limit increased from " <<
	LowEnergyLimit()/eV << " eV to " << lowEnergyLimit/eV << " eV" << G4endl;
    SetLowEnergyLimit(lowEnergyLimit);
    }

  if (HighEnergyLimit() > highEnergyLimit)
  {
    G4cout << "G4MuElecElasticModel: high energy limit decreased from " <<
        HighEnergyLimit()/MeV << " MeV to " << highEnergyLimit/MeV << " MeV" << G4endl;
    SetHighEnergyLimit(highEnergyLimit);
  }

  // Reading of data files

  G4double scaleFactor = 1e-18 * cm * cm;

  G4String fileElectron("microelec/sigma_elastic_e_Si");

  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
  G4String electron;

    // For total cross section

    electron = electronDef->GetParticleName();

    tableFile[electron] = fileElectron;

    G4MuElecCrossSectionDataSet* tableE = new G4MuElecCrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableE->LoadData(fileElectron);
    tableData[electron] = tableE;

    // For final state

    char *path = std::getenv("G4LEDATA");

    if (!path)
    {
      G4Exception("G4MuElecElasticModel::Initialise","em0006",FatalException,"G4LEDATA environment variable not set.");
      return;
    }

    std::ostringstream eFullFileName;
    eFullFileName << path << "/microelec/sigmadiff_elastic_e_Si.dat";
    std::ifstream eDiffCrossSection(eFullFileName.str().c_str());

    if (!eDiffCrossSection)
	G4Exception("G4MuElecElasticModel::Initialise","em0003",FatalException,"Missing data file: /microelec/sigmadiff_elastic_e_Si.dat");

    eTdummyVec.push_back(0.);

    while(!eDiffCrossSection.eof())
    {
	double tDummy;
	double eDummy;
	eDiffCrossSection>>tDummy>>eDummy;

	// SI : mandatory eVecm initialization
        if (tDummy != eTdummyVec.back())
        {
          eTdummyVec.push_back(tDummy);
          eVecm[tDummy].push_back(0.);
        }

        eDiffCrossSection>>eDiffCrossSectionData[tDummy][eDummy];

	// SI : only if not end of file reached !
        if (!eDiffCrossSection.eof()) eDiffCrossSectionData[tDummy][eDummy]*=scaleFactor;

        if (eDummy != eVecm[tDummy].back()) eVecm[tDummy].push_back(eDummy);

    }

    // End final state

  if (verboseLevel > 2)
    G4cout << "Loaded cross section files for MuElec Elastic model" << G4endl;

  if( verboseLevel>0 )
  {
    G4cout << "MuElec Elastic model is initialized " << G4endl
           << "Energy range: "
           << LowEnergyLimit() / eV << " eV - "
           << HighEnergyLimit() / keV << " keV"
           << G4endl;
  }

  if (isInitialised) { return; }
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;

  // InitialiseElementSelectors(particle,cuts);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuElecElasticModel::CrossSectionPerVolume(const G4Material* material,
					   const G4ParticleDefinition* p,
					   G4double ekin,
					   G4double,
					   G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4MuElecElasticModel" << G4endl;

 // Calculate total cross section for model

 G4double sigma=0;

 G4double density = material->GetTotNbOfAtomsPerVolume();

 if (material == nistSi || material->GetBaseMaterial() == nistSi)
 {
  const G4String& particleName = p->GetParticleName();

  if (ekin < highEnergyLimit)
  {
      //SI : XS must not be zero otherwise sampling of secondaries method ignored
      if (ekin < lowEnergyLimitOfModel) ekin = lowEnergyLimitOfModel;
      //

	std::map< G4String,G4MuElecCrossSectionDataSet*,std::less<G4String> >::iterator pos;
	pos = tableData.find(particleName);

	if (pos != tableData.end())
	{
	  G4MuElecCrossSectionDataSet* table = pos->second;
	  if (table != 0)
	  {
	    sigma = table->FindValue(ekin);
	  }
	}
	else
	{
	    G4Exception("G4MuElecElasticModel::ComputeCrossSectionPerVolume","em0002",FatalException,"Model not applicable to particle type.");
	}
  }

  if (verboseLevel > 3)
  {
    G4cout << "---> Kinetic energy(eV)=" << ekin/eV << G4endl;
    G4cout << " - Cross section per Si atom (cm^2)=" << sigma/cm/cm << G4endl;
    G4cout << " - Cross section per Si atom (cm^-1)=" << sigma*density/(1./cm) << G4endl;
  }

 }

 return sigma*density;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuElecElasticModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
					      const G4MaterialCutsCouple* /*couple*/,
					      const G4DynamicParticle* aDynamicElectron,
					      G4double,
					      G4double)
{

  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4MuElecElasticModel" << G4endl;

  G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();

  if (electronEnergy0 < killBelowEnergy)
  {
    fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(electronEnergy0);
    return ;
  }

  if (electronEnergy0>= killBelowEnergy && electronEnergy0 < highEnergyLimit)
  {
    G4double cosTheta = RandomizeCosTheta(electronEnergy0);

    G4double phi = 2. * pi * G4UniformRand();

    G4ThreeVector zVers = aDynamicElectron->GetMomentumDirection();
    G4ThreeVector xVers = zVers.orthogonal();
    G4ThreeVector yVers = zVers.cross(xVers);

    G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
    G4double yDir = xDir;
    xDir *= std::cos(phi);
    yDir *= std::sin(phi);

    G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers));

    fParticleChangeForGamma->ProposeMomentumDirection(zPrimeVers.unit()) ;

    fParticleChangeForGamma->SetProposedKineticEnergy(electronEnergy0);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuElecElasticModel::Theta
  (G4ParticleDefinition * particleDefinition, G4double k, G4double integrDiff)
{

  G4double theta = 0.;
  G4double valueT1 = 0;
  G4double valueT2 = 0;
  G4double valueE21 = 0;
  G4double valueE22 = 0;
  G4double valueE12 = 0;
  G4double valueE11 = 0;
  G4double xs11 = 0;
  G4double xs12 = 0;
  G4double xs21 = 0;
  G4double xs22 = 0;


  if (particleDefinition == G4Electron::ElectronDefinition())
  {
    std::vector<double>::iterator t2 = std::upper_bound(eTdummyVec.begin(),eTdummyVec.end(), k);
    std::vector<double>::iterator t1 = t2-1;

    std::vector<double>::iterator e12 = std::upper_bound(eVecm[(*t1)].begin(),eVecm[(*t1)].end(), integrDiff);
    std::vector<double>::iterator e11 = e12-1;

    std::vector<double>::iterator e22 = std::upper_bound(eVecm[(*t2)].begin(),eVecm[(*t2)].end(), integrDiff);
    std::vector<double>::iterator e21 = e22-1;

    valueT1  =*t1;
    valueT2  =*t2;
    valueE21 =*e21;
    valueE22 =*e22;
    valueE12 =*e12;
    valueE11 =*e11;

    xs11 = eDiffCrossSectionData[valueT1][valueE11];
    xs12 = eDiffCrossSectionData[valueT1][valueE12];
    xs21 = eDiffCrossSectionData[valueT2][valueE21];
    xs22 = eDiffCrossSectionData[valueT2][valueE22];

}

  if (xs11==0 || xs12==0 ||xs21==0 ||xs22==0) return (0.);

  theta = QuadInterpolator(  valueE11, valueE12,
    			     valueE21, valueE22,
			     xs11, xs12,
			     xs21, xs22,
			     valueT1, valueT2,
			     k, integrDiff );

  return theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuElecElasticModel::LinLogInterpolate(G4double e1,
						        G4double e2,
						        G4double e,
						        G4double xs1,
						        G4double xs2)
{
  G4double d1 = std::log(xs1);
  G4double d2 = std::log(xs2);
  G4double value = G4Exp(d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuElecElasticModel::LogLogInterpolate(G4double e1,
						        G4double e2,
						        G4double e,
						        G4double xs1,
						        G4double xs2)
{
  G4double a = (std::log10(xs2)-std::log10(xs1)) / (std::log10(e2)-std::log10(e1));
  G4double b = std::log10(xs2) - a*std::log10(e2);
  G4double sigma = a*std::log10(e) + b;
  G4double value = (std::pow(10.,sigma));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuElecElasticModel::QuadInterpolator(G4double e11, G4double e12,
						       G4double e21, G4double e22,
						       G4double xs11, G4double xs12,
						       G4double xs21, G4double xs22,
						       G4double t1, G4double t2,
						       G4double t, G4double e)
{
// Lin-Log
  G4double interpolatedvalue1 = LinLogInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LinLogInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LinLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuElecElasticModel::RandomizeCosTheta(G4double k)
{
  G4double integrdiff=0;
 G4double uniformRand=G4UniformRand();
 integrdiff = uniformRand;

 G4double theta=0.;
 G4double cosTheta=0.;
 theta = Theta(G4Electron::ElectronDefinition(),k/eV,integrdiff);

 cosTheta= std::cos(theta*pi/180);

 return cosTheta;
}
