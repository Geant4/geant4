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
// Created on 2016/01/18
//
// Authors: D. Sakata, W.G. Shin, S. Incerti
//
// Based on a recent release of the ELSEPA code 
// developed and provided kindly by F. Salvat et al. 
// See
// Computer Physics Communications, 165(2), 157-190. (2005)
// http://dx.doi.org/10.1016/j.cpc.2004.09.006
//

#include "G4DNAELSEPAElasticModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAMolecularMaterial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAELSEPAElasticModel::G4DNAELSEPAElasticModel(const G4ParticleDefinition*,
const G4String& nam) :
G4VEmModel(nam), isInitialised(false)
{
  verboseLevel = 0;

  G4ProductionCutsTable* theCoupleTable =
  G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = theCoupleTable->GetTableSize();

  for(G4int i=0; i<numOfCouples; ++i)
  {
    const G4MaterialCutsCouple* couple =
         theCoupleTable->GetMaterialCutsCouple(i);
    const G4Material* material = couple->GetMaterial();
    G4int nelm = material->GetNumberOfElements();
    const G4ElementVector* theElementVector = material->GetElementVector();

    if(nelm==1)
    {// Protection: only for single element
      G4int Z = 79;
      Z =  G4lrint((*theElementVector)[0]->GetZ());
      // Protection: only for GOLD
      if (Z==79){
        killBelowEnergy = 10. * eV;  // Kills e- tracking
        lowEnergyLimit  = 0   * eV;  // Must stay at zero for killing
        highEnergyLimit = 1   * GeV; // Default
        SetLowEnergyLimit (lowEnergyLimit);
        SetHighEnergyLimit(highEnergyLimit);
      }else{
        //continue;
      }
    }else{// Protection: H2O only is available
      if(material->GetName()=="G4_WATER"){
        lowEnergyLimit  = 10. * eV;  
        highEnergyLimit = 1   * MeV; 
        SetLowEnergyLimit (lowEnergyLimit);
        SetHighEnergyLimit(highEnergyLimit);
      }else{
        //continue;
      }
    }

    if (verboseLevel > 0)
    {
      G4cout << "ELSEPA Elastic model is constructed for " 
      << material->GetName() << G4endl 
      << "Energy range: "
      << lowEnergyLimit / eV << " eV - "
      << highEnergyLimit / MeV << " MeV"
      << G4endl;
    }
  }


  fParticleChangeForGamma = 0;
  fpMolDensity = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAELSEPAElasticModel::~G4DNAELSEPAElasticModel()
{
  std::map<G4int,G4DNACrossSectionDataSet*,
           std::less<G4String>>::iterator posZ;
  for (posZ = tableZData.begin(); posZ != tableZData.end(); ++posZ)
  {
    G4DNACrossSectionDataSet* table = posZ->second;
    delete table;
  }

  eEdummyVecZ.clear();
  eCumZ.clear();
  fAngleDataZ.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAELSEPAElasticModel::Initialise(const G4ParticleDefinition* particle,
const G4DataVector& )
{
  if (verboseLevel > 3)
  G4cout << "Calling G4DNAELSEPAElasticModel::Initialise()" << G4endl;

  if (isInitialised) {return;}

  if(particle->GetParticleName() != "e-")
  {
    G4Exception("G4DNAELSEPAElasticModel::Initialise","em0001",
      FatalException,"Model not applicable to particle type.");
    return;
  }
 
  G4ProductionCutsTable* theCoupleTable =
  G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = theCoupleTable->GetTableSize();
  
  // UNIT OF TCS
  G4double scaleFactor = 1.*cm*cm;

  for(G4int i=0; i<numOfCouples; ++i) 
  {
    const G4MaterialCutsCouple* couple = 
         theCoupleTable->GetMaterialCutsCouple(i);
    const G4Material* material = couple->GetMaterial();
    const G4ElementVector* theElementVector = material->GetElementVector();

    G4int nelm = material->GetNumberOfElements();
    if (nelm==1){// Protection: only for single element
      G4int Z =  G4lrint((*theElementVector)[0]->GetZ());
      if (Z!=79)// Protection: only for GOLD
      {
        continue;
      }
      
      if (Z>0) 
      {
        G4String fileZElectron("dna/sigma_elastic_e_elsepa_Z");
        std::ostringstream oss;
        oss.str("");
        oss.clear(stringstream::goodbit);
        oss << Z;
        fileZElectron += oss.str()+"_muffintin";
        
        G4DNACrossSectionDataSet* tableZE =
          new G4DNACrossSectionDataSet
            (new G4LogLogInterpolation, eV,scaleFactor );
        tableZE->LoadData(fileZElectron);
        tableZData[Z] = tableZE;
      
        std::ostringstream eFullFileNameZ;
        char *path = getenv("G4LEDATA");
        if (!path)
        {
          G4Exception("G4DNAELSEPAElasticModel::Initialise","em0002",
            FatalException,"G4LEDATA environment variable not set.");
          return;
        }

        eFullFileNameZ.str("");
        eFullFileNameZ.clear(stringstream::goodbit);
      
        eFullFileNameZ 
          << path 
          << "/dna/sigmadiff_cumulated_elastic_e_elsepa_Z" 
          << Z << "_muffintin.dat";
      
        std::ifstream eDiffCrossSectionZ(eFullFileNameZ.str().c_str());
      
        if (!eDiffCrossSectionZ)
        {
          G4Exception("G4DNAELSEPAElasticModel::Initialise","em0003",
            FatalException,"Missing data file for cumulated DCS");
          return;
        }
      
        eEdummyVecZ.clear();
        eCumZ.clear();
        fAngleDataZ.clear();
        
        eEdummyVecZ[Z].push_back(0.);
        do
        {
          G4double eDummy;
          G4double cumDummy;
          eDiffCrossSectionZ>>eDummy>>cumDummy;
          if (eDummy != eEdummyVecZ[Z].back())
          {
           eEdummyVecZ[Z].push_back(eDummy);
           eCumZ[Z][eDummy].push_back(0.);
          }
          eDiffCrossSectionZ>>fAngleDataZ[Z][eDummy][cumDummy];
          if (cumDummy != eCumZ[Z][eDummy].back())
          {
            eCumZ[Z][eDummy].push_back(cumDummy);
          }
        }while(!eDiffCrossSectionZ.eof());
      } 

    }else{// Protection: H2O only is available
      if(material->GetName()=="G4_WATER"){
        if (LowEnergyLimit() < 10*eV)
        {
          G4cout<<"G4DNAELSEPAElasticModel: low energy limit increased from "
                << LowEnergyLimit()/eV << " eV to " << 10 << " eV"
                << G4endl;
          SetLowEnergyLimit(10.*eV);
        }

        if (HighEnergyLimit() > 1.*MeV)
        {
          G4cout<<"G4DNAELSEPAElasticModel: high energy limit decreased from "
                << HighEnergyLimit()/MeV << " MeV to " << 1. << " MeV"
                << G4endl;
          SetHighEnergyLimit(1.*MeV);
        }

        G4String fileZElectron("dna/sigma_elastic_e_elsepa_muffin");

        G4DNACrossSectionDataSet* tableZE =
          new G4DNACrossSectionDataSet(
             new G4LogLogInterpolation, eV,scaleFactor );
        tableZE->LoadData(fileZElectron);
        tableZData[0] = tableZE;

        std::ostringstream eFullFileNameZ;

        char *path = getenv("G4LEDATA");
        if (!path)
        {
          G4Exception("G4DNAELSEPAElasticModel::Initialise","em0004",
            FatalException,"G4LEDATA environment variable not set.");
          return;
        }

        eFullFileNameZ.str("");
        eFullFileNameZ.clear(stringstream::goodbit);

        eFullFileNameZ
          << path
          <<  "/dna/sigmadiff_cumulated_elastic_e_elsepa_muffin.dat";

        std::ifstream eDiffCrossSectionZ(eFullFileNameZ.str().c_str());

        if (!eDiffCrossSectionZ)
         G4Exception("G4DNAELSEPAElasticModel::Initialise","em0005",
         FatalException,
         "Missing data file for cumulated DCS");

        eEdummyVecZ.clear();
        eCumZ.clear();
        fAngleDataZ.clear();

        eEdummyVecZ[0].push_back(0.);

        do
        {
          G4double eDummy;
          G4double cumDummy;
          eDiffCrossSectionZ>>eDummy>>cumDummy;
          if (eDummy != eEdummyVecZ[0].back())
          {
           eEdummyVecZ[0].push_back(eDummy);
           eCumZ[0][eDummy].push_back(0.);
          }
          eDiffCrossSectionZ>>fAngleDataZ[0][eDummy][cumDummy];
          if (cumDummy != eCumZ[0][eDummy].back())
         eCumZ[0][eDummy].push_back(cumDummy);
        }while(!eDiffCrossSectionZ.eof());
      }
    }
    if (verboseLevel > 2)
    G4cout << "Loaded cross section files of ELSEPA Elastic model for"
           << material->GetName() << G4endl;

    if( verboseLevel>0 )
    {
      G4cout << "ELSEPA elastic model is initialized " << G4endl
      << "Energy range: "
      << LowEnergyLimit() /  eV << " eV - "
      << HighEnergyLimit()/ MeV << " MeV"
      << G4endl;
    }
  } // Loop on couples

  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAELSEPAElasticModel::CrossSectionPerVolume
(const G4Material* material,
 const G4ParticleDefinition* particle,
 G4double ekin,
 G4double,
 G4double)
{

  if (verboseLevel > 3)
  {
    G4cout <<
    "Calling CrossSectionPerVolume() of G4DNAELSEPAElasticModel"
    << G4endl;
  }

  G4double atomicNDensity=0.0;
  G4double sigma=0;

  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nelm = material->GetNumberOfElements();
  if (nelm==1) {// Protection: only for single element
    // Protection: only for GOLD
    if (material->GetZ()!=79) return 0;

    G4int Z = G4lrint((*theElementVector)[0]->GetZ());

    const G4String& particleName = particle->GetParticleName();
    atomicNDensity = material->GetAtomicNumDensityVector()[0];
    if(atomicNDensity!= 0.0)
    {
      if (ekin < highEnergyLimit)
      {
        if (ekin < killBelowEnergy) return DBL_MAX;

        std::map< G4int,G4DNACrossSectionDataSet*,
                  std::less<G4String> >::iterator pos;
        pos = tableZData.find(Z);
        
        if (pos != tableZData.end())
        {
          G4DNACrossSectionDataSet* table = pos->second;
          if (table != 0)
          {
            // XS takes its 10 eV value below 10 eV for GOLD
            if (ekin < 10*eV) sigma = table->FindValue(10*eV);
            else sigma = table->FindValue(ekin);
          }
        }
        else
        {
          G4Exception("G4DNAELSEPAElasticModel::ComputeCrossSectionPerVolume",
            "em0006",FatalException,"Model not applicable to particle type.");
        }
      }
    }
    if (verboseLevel > 2)
    {
      G4cout << "__________________________________" << G4endl;
      G4cout << "=== G4DNAELSEPAElasticModel - XS INFO START" << G4endl;
      G4cout << "=== Material is made of one element with Z =" << Z << G4endl;
      G4cout << "=== Kinetic energy(eV)=" << ekin/eV << " particle : " 
             << particleName << G4endl;
      G4cout << "=== Cross section per atom for Z="<<Z<<" is (cm^2)" 
             << sigma/cm/cm << G4endl;
      G4cout << "=== Cross section per atom for Z="<<Z<<" is (cm^-1)=" 
             << sigma*atomicNDensity/(1./cm) << G4endl;
      G4cout << "=== G4DNAELSEPAElasticModel - XS INFO END" << G4endl;
    }
  }else{
    fpMolDensity =
     G4DNAMolecularMaterial::Instance()->
     GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));
     atomicNDensity = (*fpMolDensity)[material->GetIndex()];
     if(atomicNDensity!= 0.0)
     {
       if (ekin < HighEnergyLimit() && ekin >= LowEnergyLimit())
       {
         std::map< G4int,G4DNACrossSectionDataSet*,
         std::less<G4String> >::iterator pos;
         pos = tableZData.find(0); // the data is stored as Z=0
         //SI : XS must not be zero 
        //     otherwise sampling of secondaries method ignored
         if (pos != tableZData.end())
         {
           G4DNACrossSectionDataSet* table = pos->second;
           if (table != 0)
           {
             sigma = table->FindValue(ekin);
           }
         }
       }
     }
     if (verboseLevel > 2)
     {
       G4cout << "__________________________________" << G4endl;
       G4cout << "=== G4DNAELSEPAElasticModel - XS INFO START" << G4endl;
       G4cout << "=== Kinetic energy(eV)=" << ekin/eV 
              << " particle : " << particle->GetParticleName() << G4endl;
       G4cout << "=== Cross section per water molecule (cm^2)=" 
              << sigma/cm/cm << G4endl;
       G4cout << "=== Cross section per water molecule (cm^-1)=" 
              << sigma*atomicNDensity/(1./cm) << G4endl;
       G4cout << "=== G4DNAELSEPAElasticModel - XS INFO END" << G4endl;
     }
  }

  return sigma*atomicNDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAELSEPAElasticModel::SampleSecondaries(
      std::vector<G4DynamicParticle*>*,
      const G4MaterialCutsCouple* couple,
      const G4DynamicParticle* aDynamicElectron,
      G4double,
      G4double)
{

  if (verboseLevel > 3){
    G4cout << 
    "Calling SampleSecondaries() of G4DNAELSEPAElasticModel" 
    << G4endl;
  }

  G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();

  const G4Material* material = couple->GetMaterial();
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nelm = material->GetNumberOfElements();
  if (nelm==1){// Protection: only for single element
    G4int Z =  G4lrint((*theElementVector)[0]->GetZ());
    if (Z!=79) return;
    if (electronEnergy0 < killBelowEnergy)
    {
      fParticleChangeForGamma->SetProposedKineticEnergy(0.);
      fParticleChangeForGamma->ProposeMomentumDirection(G4ThreeVector(0,0,0));
      fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(electronEnergy0);
      return;
    }

    if(electronEnergy0>= killBelowEnergy && electronEnergy0 < highEnergyLimit)
    {
      G4double cosTheta = 0;
      if (electronEnergy0>=10*eV){
         cosTheta = RandomizeCosTheta(Z,electronEnergy0);
      }else { 
        cosTheta = RandomizeCosTheta(Z,10*eV);
      }

      G4double phi = 2. * CLHEP::pi * G4UniformRand();
      
      G4ThreeVector zVers = aDynamicElectron->GetMomentumDirection();
      G4ThreeVector xVers = zVers.orthogonal();
      G4ThreeVector yVers = zVers.cross(xVers);

      G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
      G4double yDir = xDir;
      xDir *= std::cos(phi);
      yDir *= std::sin(phi);

      G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers));
      fParticleChangeForGamma->ProposeMomentumDirection(zPrimeVers.unit());
      fParticleChangeForGamma->SetProposedKineticEnergy(electronEnergy0);
      
    }
  }else{
    if(material->GetName()=="G4_WATER"){
      //The data for water is stored as Z=0
      G4double cosTheta = RandomizeCosTheta(0,electronEnergy0);

      G4double phi = 2. * pi * G4UniformRand();

      G4ThreeVector zVers = aDynamicElectron->GetMomentumDirection();
      G4ThreeVector xVers = zVers.orthogonal();
      G4ThreeVector yVers = zVers.cross(xVers);

      G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
      G4double yDir = xDir;
      xDir *= std::cos(phi);
      yDir *= std::sin(phi);

      G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers));
      fParticleChangeForGamma->ProposeMomentumDirection(zPrimeVers.unit());
      fParticleChangeForGamma->SetProposedKineticEnergy(electronEnergy0);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAELSEPAElasticModel::Theta(G4int Z,
         G4ParticleDefinition * particleDefinition,
         G4double k,
         G4double integrDiff)
{

 G4double theta   = 0.;
 G4double valueE1 = 0.;
 G4double valueE2 = 0.;
 G4double valuecum21 = 0.;
 G4double valuecum22 = 0.;
 G4double valuecum12 = 0.;
 G4double valuecum11 = 0.;
 G4double a11 = 0.;
 G4double a12 = 0.;
 G4double a21 = 0.;
 G4double a22 = 0.;

 if (particleDefinition == G4Electron::ElectronDefinition())
 {
  std::vector<G4double>::iterator e2 
             = std::upper_bound(eEdummyVecZ[Z].begin(),
             eEdummyVecZ[Z].end(), k);
  std::vector<G4double>::iterator e1 = e2 - 1;

  std::vector<G4double>::iterator cum12 
             = std::upper_bound(eCumZ[Z][(*e1)].begin(),
             eCumZ[Z][(*e1)].end(),integrDiff);
  
  std::vector<G4double>::iterator cum11 = cum12 - 1;

  std::vector<G4double>::iterator cum22 
             = std::upper_bound(eCumZ[Z][(*e2)].begin(),
             eCumZ[Z][(*e2)].end(),integrDiff);
  
  std::vector<G4double>::iterator cum21 = cum22 - 1;

  valueE1  = *e1;
  valueE2  = *e2;
  valuecum11 = *cum11;
  valuecum12 = *cum12;
  valuecum21 = *cum21;
  valuecum22 = *cum22;


  a11 = fAngleDataZ[Z][valueE1][valuecum11];
  a12 = fAngleDataZ[Z][valueE1][valuecum12];
  a21 = fAngleDataZ[Z][valueE2][valuecum21];
  a22 = fAngleDataZ[Z][valueE2][valuecum22];

 }

 if (a11 == 0 && a12 == 0 && a21 == 0 && a22 == 0) return (0.);

 theta = QuadInterpolator(valuecum11, valuecum12, valuecum21, valuecum22, 
          a11, a12,a21, a22, valueE1, valueE2, k, integrDiff);
 return theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//
G4double G4DNAELSEPAElasticModel::LogLinInterpolate(G4double e1,
G4double e2,
G4double e,
G4double xs1,
G4double xs2)
{
 G4double value=0.;
 if(e1!=0){
   G4double a = std::log10(e)  - std::log10(e1);
   G4double b = std::log10(e2) - std::log10(e);
   value = xs1 + a/(a+b)*(xs2-xs1);
 }
 else{
   G4double d1 = xs1;
   G4double d2 = xs2;
   value = (d1 + (d2 - d1) * (e - e1) / (e2 - e1));
 }

 return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAELSEPAElasticModel::LinLogInterpolate(G4double e1,
G4double e2,
G4double e,
G4double xs1,
G4double xs2)
{
 G4double d1 = std::log10(xs1);
 G4double d2 = std::log10(xs2);
 G4double value = std::pow(10,(d1 + (d2 - d1) * (e - e1) / (e2 - e1)));
 return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAELSEPAElasticModel::LinLinInterpolate(G4double e1,
G4double e2,
G4double e,
G4double xs1,
G4double xs2)
{
 G4double d1 = xs1;
 G4double d2 = xs2;
 G4double value = (d1 + (d2 - d1) * (e - e1) / (e2 - e1));
 return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAELSEPAElasticModel::LogLogInterpolate(G4double e1,
G4double e2,
G4double e,
G4double xs1,
G4double xs2)
{
 G4double a = (std::log10(xs2) - std::log10(xs1))
             / (std::log10(e2) - std::log10(e1));
 G4double b = std::log10(xs2) - a * std::log10(e2);
 G4double sigma = a * std::log10(e) + b;
 G4double value = (std::pow(10., sigma));
 return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAELSEPAElasticModel::QuadInterpolator(
G4double cum11,
G4double cum12,
G4double cum21,
G4double cum22,
G4double a11,
G4double a12,
G4double a21,
G4double a22,
G4double e1,
G4double e2,
G4double t,
G4double cum)
{
   G4double value=0;
   G4double interpolatedvalue1=0;
   G4double interpolatedvalue2=0;

   if(cum11!=0){
     interpolatedvalue1 = LinLogInterpolate(cum11, cum12, cum, a11, a12);
   }
   else{
     interpolatedvalue1 = LinLinInterpolate(cum11, cum12, cum, a11, a12);
   }
   if(cum21!=0){
     interpolatedvalue2 = LinLogInterpolate(cum21, cum22, cum, a21, a22);
   }
   else{
     interpolatedvalue2 = LinLinInterpolate(cum21, cum22, cum, a21, a22);
   }

   value = LogLinInterpolate(e1,e2,t,interpolatedvalue1,interpolatedvalue2);

 return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAELSEPAElasticModel::RandomizeCosTheta(G4int Z, G4double k)
{

  G4double integrdiff = 0.;
  G4double uniformRand = G4UniformRand();
  integrdiff = uniformRand;

  G4double theta = 0.;
  G4double cosTheta = 0.;
  theta = Theta(Z, G4Electron::ElectronDefinition(), k / eV, integrdiff);

  cosTheta = std::cos(theta * CLHEP::pi / 180.); 

  return cosTheta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAELSEPAElasticModel::SetKillBelowThreshold(G4double threshold)
{
  killBelowEnergy = threshold;

  if (threshold < 10 * eV)
  {
    G4cout<< "*** WARNING : the G4DNAELSEPAElasticModel model is not "
    "defined below 10 eV !" << G4endl;
  }
}
