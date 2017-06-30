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
// CPA100 elastic model class for electrons
//
// Based on the work of M. Terrissol and M. C. Bordage
//
// Users are requested to cite the following papers:
// - M. Terrissol, A. Baudre, Radiat. Prot. Dosim. 31 (1990) 175-177
// - M.C. Bordage, J. Bordes, S. Edel, M. Terrissol, X. Franceries, 
//   M. Bardies, N. Lampe, S. Incerti, Phys. Med. 32 (2016) 1833-1840
//
// Authors of this class: 
// M.C. Bordage, M. Terrissol, S. Edel, J. Bordes, S. Incerti
//
// 15.01.2014: creation
//

#include "G4DNACPA100ElasticModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAMolecularMaterial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

// #define CPA100_VERBOSE

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACPA100ElasticModel::G4DNACPA100ElasticModel(const G4ParticleDefinition*,
                                             const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{

  //killBelowEnergy = 11. * eV; // Default
  //killBelowEnergy = 10.481 * eV; 
  //lowEnergyLimit = 11 * eV; 
  //highEnergyLimit = 255955 * eV;
  SetLowEnergyLimit(11*eV);
  SetHighEnergyLimit(255955 * eV);

  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

#ifdef UEHARA_VERBOSE  
  if( verboseLevel>0 ) 
  { 
    G4cout << "CPA100 Elastic model is constructed " << G4endl
           << "Energy range: "
           << lowEnergyLimit / eV << " eV - "
           << highEnergyLimit / keV << " keV"
           << G4endl;
  }
#endif
  
  fParticleChangeForGamma = 0;
  fpMolWaterDensity = 0;

  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACPA100ElasticModel::~G4DNACPA100ElasticModel()
{  
  // For total cross section
  
  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4DNACrossSectionDataSet* table = pos->second;
    delete table;
  }

  // For final state
   
  eVecm.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNACPA100ElasticModel::Initialise(const G4ParticleDefinition* 
                                       particle,
                                       const G4DataVector& /*cuts*/)
{

#ifdef UEHARA_VERBOSE
  if (verboseLevel > 3)
    G4cout << "Calling G4DNACPA100ElasticModel::Initialise()" << G4endl;
#endif

  if(particle->GetParticleName() != "e-")
  {
    G4Exception("*** WARNING: the G4DNACPA100ElasticModel is "
                "not intented to be used with another particle than the electron",
                "",FatalException,"") ;
  }
  
  // Energy limits
  
  if (LowEnergyLimit() < 11.*eV)
  {
    G4cout << "G4DNACPA100ElasticModel: low energy limit increased from " << 
    LowEnergyLimit()/eV << " eV to " << 11 << " eV" << G4endl;
    SetLowEnergyLimit(11.*eV);
  }

  if (HighEnergyLimit() > 255955.*eV)
  {
    G4cout << "G4DNACPA100ElasticModel: high energy limit decreased from " << 
        HighEnergyLimit()/keV << " keV to " << 255.955 << " keV" 
        << G4endl;
    SetHighEnergyLimit(255955.*eV);
  }

  // Reading of data files 
  
  G4double scaleFactor = 1e-20*m*m;

  G4String fileElectron("dna/sigma_elastic_e_cpa100");

  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
  G4String electron;
 
  // *** ELECTRON
  
    // For total cross section
    
    electron = electronDef->GetParticleName();

    tableFile[electron] = fileElectron;

    G4DNACrossSectionDataSet* tableE = 
     new G4DNACrossSectionDataSet(new G4LogLogInterpolation, 
     eV,scaleFactor );
     
    /*
    G4DNACrossSectionDataSet* tableE = 
      new G4DNACrossSectionDataSet(new G4DNACPA100LogLogInterpolation, 
            eV,scaleFactor ); 
    */
    
    tableE->LoadData(fileElectron);
    
    tableData[electron] = tableE;
    
    // For final state

    char *path = getenv("G4LEDATA");
 
    if (!path)
    {
      G4Exception("G4DNACPA100ElasticModel::Initialise","em0006",
                  FatalException,"G4LEDATA environment variable not set.");
      return;
    }
       
    std::ostringstream eFullFileName;

    eFullFileName << path 
    << "/dna/sigmadiff_cumulated_elastic_e_cpa100.dat";
    std::ifstream eDiffCrossSection(eFullFileName.str().c_str());
     
    if (!eDiffCrossSection) 
      G4Exception("G4DNACPA100ElasticModel::Initialise","em0003",
                  FatalException,
      "Missing data file:/dna/sigmadiff_cumulated_elastic_e_cpa100.dat");
    
    // March 25th, 2014 - Vaclav Stepan, Sebastien Incerti
    // Added clear for MT

    eTdummyVec.clear();
    eVecm.clear();
    eDiffCrossSectionData.clear();

    //

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

        if (eDummy != eVecm[tDummy].back()) eVecm[tDummy].push_back(eDummy);
 
    }

    // End final state
    
#ifdef UEHARA_VERBOSE
  if (verboseLevel > 2) 
    G4cout << "Loaded cross section files for CPA100 Elastic model" << G4endl;
#endif

#ifdef UEHARA_VERBOSE
  if( verboseLevel>0 ) 
  { 
    G4cout << "CPA100 Elastic model is initialized " << G4endl
           << "Energy range: "
           << LowEnergyLimit() / eV << " eV - "
           << HighEnergyLimit() / keV << " keV"
           << G4endl;
  }
#endif

  // Initialize water density pointer
  fpMolWaterDensity = G4DNAMolecularMaterial::Instance()
    ->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  if (isInitialised) { return; }
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ElasticModel::CrossSectionPerVolume
(const G4Material* material,
        const G4ParticleDefinition* p,
        G4double ekin,
        G4double,
        G4double)
{

#ifdef UEHARA_VERBOSE
 if (verboseLevel > 3)
    G4cout << 
    "Calling CrossSectionPerVolume() of G4DNACPA100ElasticModel" << G4endl;
#endif

 // Calculate total cross section for model

 G4double sigma=0;

 G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

 if(waterDensity!= 0.0)
  {
  const G4String& particleName = p->GetParticleName();

  if (ekin < HighEnergyLimit() && ekin >= LowEnergyLimit())
  {
      //SI : XS must not be zero otherwise sampling of secondaries 
      // method ignored
      
      std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
      pos = tableData.find(particleName);
 
      if (pos != tableData.end())
      {
        G4DNACrossSectionDataSet* table = pos->second;
        if (table != 0)
        {
          sigma = table->FindValue(ekin);

 //
 //Dump in non-MT mode
 //
 /*
        G4double minEnergy = 10.481  * eV;
        G4double maxEnergy = 255955. * eV;
        G4int nEnergySteps = 1000;
        G4double energy(minEnergy);
        G4double 
          stpEnergy(std::pow(maxEnergy/energy, 
          1./static_cast<G4double>(nEnergySteps-1)));
        G4int step(nEnergySteps);
        system ("rm -rf elastic-cpa100.out");
        FILE* myFile=fopen("elastic-cpa100.out","a");
        while (step>0)
        {
          step--;
          fprintf (myFile,"%16.9le %16.9le\n",
            energy/eV,
            table->FindValue(energy)/(1e-20*m*m));
          energy*=stpEnergy;
        }
        fclose (myFile);
        abort();
 */
 //
 // end of dump
 //
        
              
   }
 }
 else
 {
    G4Exception("G4DNACPA100ElasticModel::ComputeCrossSectionPerVolume",
      "em0002",
      FatalException,"Model not applicable to particle type.");
 }
  
  }

#ifdef UEHARA_VERBOSE
  if (verboseLevel > 2)
  {
    G4cout << "__________________________________" << G4endl;
    G4cout << "G4DNACPA100ElasticModel - XS INFO START" << G4endl;
    G4cout << "Kinetic energy(eV)=" << ekin/eV << " particle : " << particleName << G4endl;
    G4cout << "Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
    G4cout << "Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
    // G4cout << " - Cross section per water molecule (cm^-1)=" 
    // << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
    G4cout << "G4DNACPA100ElasticModel - XS INFO END" << G4endl;
  } 
#endif

 } 
  
   return sigma*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNACPA100ElasticModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
           const G4MaterialCutsCouple* /*couple*/,
           const G4DynamicParticle* aDynamicElectron,
           G4double,
           G4double)
{
#ifdef UEHARA_VERBOSE
  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4DNACPA100ElasticModel" << G4endl;
#endif

  G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();
  
  {  
    G4double cosTheta = RandomizeCosTheta(electronEnergy0);
    G4double phi = 2. * pi * G4UniformRand();

    G4ThreeVector zVers = aDynamicElectron->GetMomentumDirection();

    //G4ThreeVector xVers = zVers.orthogonal();
    //G4ThreeVector yVers = zVers.cross(xVers);
    //G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
    //G4double yDir = xDir;
    //xDir *= std::cos(phi);
    //yDir *= std::sin(phi);

    // Computation of scattering angles (from Subroutine DIRAN in CPA100)

    G4double CT1, ST1, CF1, SF1, CT2, ST2, CF2, SF2;
    G4double sinTheta = std::sqrt (1-cosTheta*cosTheta);

    CT1=0;
    ST1=0;
    CF1=0;
    SF1=0;
    CT2=0;
    ST2=0;
    CF2=0;
    SF2=0;

    CT1 = zVers.z();
    ST1=std::sqrt(1.-CT1*CT1);

    if (ST1!=0) CF1 = zVers.x()/ST1; else CF1 = std::cos(2. * pi * G4UniformRand());
    if (ST1!=0) SF1 = zVers.y()/ST1; else SF1 = std::sqrt(1.-CF1*CF1);

    G4double A3, A4, A5, A2, A1;
    A3=0;
    A4=0;
    A5=0;
    A2=0;
    A1=0;

    A3 = sinTheta*std::cos(phi);
    A4 = A3*CT1 + ST1*cosTheta;
    A5 = sinTheta * std::sin(phi);
    A2 = A4 * SF1 + A5 * CF1;
    A1 = A4 * CF1 - A5 * SF1;

    CT2 = CT1*cosTheta - ST1*A3;
    ST2 = std::sqrt(1.-CT2*CT2);

    if (ST2==0) ST2=1E-6;
    CF2 = A1/ST2;
    SF2 = A2/ST2;

    /*
    G4cout << "CT1=" << CT1 << G4endl;
    G4cout << "ST1=" << ST1 << G4endl;
    G4cout << "CF1=" << CF1 << G4endl;
    G4cout << "SF1=" << SF1 << G4endl;
    G4cout << "cosTheta=" << cosTheta << G4endl;
    G4cout << "sinTheta=" << sinTheta << G4endl;
    G4cout << "cosPhi=" << std::cos(phi) << G4endl;
    G4cout << "sinPhi=" << std::sin(phi) << G4endl;
    G4cout << "CT2=" << CT2 << G4endl;
    G4cout << "ST2=" << ST2 << G4endl;
    G4cout << "CF2=" << CF2 << G4endl;
    G4cout << "SF2=" << SF2 << G4endl;
    */

    G4ThreeVector zPrimeVers(ST2*CF2,ST2*SF2,CT2);

    //

    fParticleChangeForGamma->ProposeMomentumDirection(zPrimeVers.unit()) ;

    //fParticleChangeForGamma->SetProposedKineticEnergy(electronEnergy0);
    
    if (!statCode) 
    
    fParticleChangeForGamma->SetProposedKineticEnergy
                      (electronEnergy0-1.214E-4*(1.-cosTheta)*electronEnergy0);
            
    else fParticleChangeForGamma->SetProposedKineticEnergy(electronEnergy0);

    fParticleChangeForGamma->ProposeLocalEnergyDeposit(1.214E-4*(1.-cosTheta)*electronEnergy0);
        
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ElasticModel::Theta
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
 
    std::vector<double>::iterator e12 = std::upper_bound(eVecm[(*t1)].begin(),eVecm[(*t1)].end(), 
     integrDiff);
    std::vector<double>::iterator e11 = e12-1;
   
    std::vector<double>::iterator e22 = std::upper_bound(eVecm[(*t2)].begin(),eVecm[(*t2)].end(),
     integrDiff);    
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
    
    //FOR CPA100
    //if(k==valueT1) xs22 = eDiffCrossSectionData[valueT1][valueE12];

  }
  
  if (xs11==0 && xs12==0 && xs21==0 && xs22==0) return (0.);   

  // FOR CPA100
  
  
  theta = QuadInterpolator   ( valueE11, valueE12, 
              valueE21, valueE22, 
          xs11, xs12, 
          xs21, xs22, 
          valueT1, valueT2, 
            k, integrDiff );
  
  return theta;
  

  //FOR CPA100  
  //return xs22;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ElasticModel::LinLogInterpolate(G4double e1, 
              G4double e2, 
              G4double e, 
              G4double xs1, 
              G4double xs2)
{
  G4double d1 = std::log(xs1);
  G4double d2 = std::log(xs2);
  G4double value = std::exp(d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ElasticModel::LinLinInterpolate(G4double e1, 
              G4double e2, 
              G4double e, 
              G4double xs1, 
              G4double xs2)
{
  G4double d1 = xs1;
  G4double d2 = xs2;
  G4double value = (d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ElasticModel::LogLogInterpolate(G4double e1, 
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


G4double G4DNACPA100ElasticModel::QuadInterpolator(G4double e11, G4double e12, 
             G4double e21, G4double e22, 
             G4double xs11, G4double xs12, 
             G4double xs21, G4double xs22, 
             G4double t1, G4double t2, 
             G4double t, G4double e)
{
  // Log-Log
/*
  G4double interpolatedvalue1 = LogLogInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LogLogInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LogLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);


  // Lin-Log
  G4double interpolatedvalue1 = LinLogInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LinLogInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LinLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
*/

  // Lin-Lin
  G4double interpolatedvalue1 = LinLinInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LinLinInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LinLinInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100ElasticModel::RandomizeCosTheta(G4double k) 
{
 
 G4double integrdiff=0; // PROBABILITY between 0 and 1.
 G4double uniformRand=G4UniformRand();
 integrdiff = uniformRand;
 
 G4double cosTheta=0.;

 // 1 - COS THETA is read from the data file 
 cosTheta = 1 - Theta(G4Electron::ElectronDefinition(),k/eV,integrdiff);

 //
 //
 //Dump
 //
 //G4cout << "theta=" << theta << G4endl;
 //G4cout << "cos theta=" << std::cos(theta*pi/180) << G4endl;
 //G4cout << "sin theta=" << std::sin(theta*pi/180) << G4endl;
 //G4cout << "acos(cos theta)=" << std::acos(cosTheta) << G4endl;
 //G4cout << "cos theta="<< cosTheta << G4endl;
 //G4cout << "1 - cos theta="<< 1. - cosTheta << G4endl;
 //G4cout << "sin theta=" << std::sqrt(1-cosTheta*cosTheta) << G4endl;
 //
 /*
 G4double minProb = 0; // we scan probability between 0 and one
 G4double maxProb = 1;
 G4int nProbSteps = 100;
 G4double prob(minProb);
 G4double stepProb((maxProb-minProb)/static_cast<G4double>(nProbSteps));
 G4int step(nProbSteps);
 system ("rm -rf elastic-cumul-cpa100-100keV.out");
 FILE* myFile=fopen("elastic-cumul-cpa100-100keV.out","a");
 while (step>=0)
 {
    step--;
    fprintf (myFile,"%16.9le %16.9le\n",
    prob,
    Theta(G4Electron::ElectronDefinition(),100000,prob)); // SELECT NRJ IN eV !!!
    prob=prob+stepProb;
 }
 fclose (myFile);
 abort();
 */
 //
 // end of dump
 //

 return cosTheta; 
}
