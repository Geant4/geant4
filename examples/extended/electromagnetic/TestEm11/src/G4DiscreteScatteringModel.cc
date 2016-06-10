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
#include "G4DiscreteScatteringModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4AtomicShells.hh"
#include "G4Track.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ElementData.hh"
#include "G4Exp.hh"
#include "G4Log.hh"

#include <iostream>
#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;
using namespace CLHEP;

G4ElementData*  G4DiscreteScatteringModel::fCdf = nullptr;
G4ElementData*  G4DiscreteScatteringModel::fTcs = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DiscreteScatteringModel::G4DiscreteScatteringModel(G4int iNumAngles)
  : G4VEmModel("DiscrScat"), fParticleChange(nullptr),fAnalogModel("pwe"), 
    fNumAngles(iNumAngles), fLowEnergyLimit(2*keV)
{        
  SetHighEnergyLimit(100.*MeV); 
  SetLowEnergyLimit(fLowEnergyLimit); 
  if(IsMaster() && fCdf == nullptr) {
    fCdf = new G4ElementData();
    fTcs = new G4ElementData();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DiscreteScatteringModel::~G4DiscreteScatteringModel()
{
  if(IsMaster()) {
    delete fCdf;
    delete fTcs;
    fCdf = nullptr;
    fTcs = nullptr;     
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DiscreteScatteringModel::Initialise(const G4ParticleDefinition*, 
                                           const G4DataVector&)
{
  if(nullptr == fParticleChange) {
    fParticleChange = GetParticleChangeForGamma();
  }
  if(!IsMaster()) { return; }

  G4cout << "G4DiscreteScatteringModel::Initialise start"<<G4endl;

  const G4int maxZ = 100;
    
  char *path = getenv("G4GBFPDATA");
  if (!path)
    {
      G4Exception("G4DiscreteScatteringModel::Initialise","em0006",
                  FatalException,"G4GBFPDATA environment variable not set.");
      return;
    }
       
  std::ostringstream eFullFileName;
  eFullFileName << path;
  
  G4ProductionCutsTable* theCoupleTable =
    G4ProductionCutsTable::GetProductionCutsTable();

  G4int numOfCouples = theCoupleTable->GetTableSize();

  for(G4int i=0; i<numOfCouples; ++i)
    { 
      const G4Material* material = 
        theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
      const G4ElementVector* theElementVector = material->GetElementVector();
      G4int nelm = material->GetNumberOfElements();

      for (G4int j=0; j<nelm; ++j)
        {   
          G4int Z = G4lrint((*theElementVector)[j]->GetZ());
          if(Z < 1)          { Z = 1; }
          else if(Z > maxZ)  { Z = maxZ; }
          if(!fTcs->GetElementData(Z)) { ReadData(Z, path); }
        }
    }
  G4cout << "G4DiscreteScatteringModel::Initialise completed"<<G4endl;        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DiscreteScatteringModel::ReadData(G4int Z, const G4String& argFileName)
{

  // Set strings with paths to tcs data file. Currently not a variable but this
  // this will change
  G4String fullFileName(argFileName);

  stringstream ss1;//create a stringstream
  ss1 << Z;//add number to the stream
  stringstream ss2;//create a stringstream
  ss2 << fNumAngles;
  G4String tcsPath = fullFileName+"/log_gbfp_"+fAnalogModel
    +"_tcs_"+ss1.str()+"_"+ss2.str()+".dat";

  // Do some error checking and exiting if data does not exist
  std::ifstream in(tcsPath, std::ifstream::binary|std::ifstream::in);
  if (!in.is_open())
    {
      G4String message("Data file \"");
      message+=tcsPath;
      message+="\" not found";
      G4Exception("G4DiscreteScatteringModel::LoadData","em0003",
                      FatalException,message);
      return;
    }
  //G4cout<<"Reading in data file "<< tcsPath<<G4endl;

  // Create a temporary G4PhysicsVector object pointer 
  G4PhysicsVector* tempData = new G4PhysicsVector(false);
  
  // Use retrieve to read in the total cross section (tcs) data
  tempData->Retrieve(in,true);
  
  // Convert tcs from cm^2 to mm^2
  //tempData->ScaleVector(1.0, 100.0);
  
  // store pass this data to the tcs object and initialise for current element 
  fTcs->InitialiseForElement(Z,tempData);  
  in.close();
                        
  // Set strings with paths to cdf data files. Currently not a variable 
  // but this this will change
  G4String cdfPath = fullFileName+"/gbfp_"+fAnalogModel+"_cdf_"
    +ss1.str()+"_"+ss2.str()+".dat";

  // Do some error checking and exiting if data does not exist
  std::ifstream in2(cdfPath, std::ifstream::binary|std::ifstream::in);
  if (!in2.is_open())
    {
      G4String message("Data file \"");
      message+=cdfPath;
      message+="\" not found";
      G4Exception("G4DiscreteScatteringModel::LoadData","em0003",
                      FatalException,message);
      return;
    }
    
  //G4cout<<"Reading in data file "<< cdfPath<<G4endl;

  // The cumulative distribution functions (CDF) for energy E_j and X_i 
  // on (0,1) is C(E_j,X_i) and read-in/stored at this time. 
  // For the purposes of this model, each energy grid point where the CDF 
  // is evaluated is considered a component. The number of energy grid 
  // points is consistent with the fTcs data, so the following int is set 
  // by calling G4PhysicsVector::GetVectorLength().
  G4int numEnergies = fTcs->GetElementData(Z)->GetVectorLength();
  
  // The ElementData object pointer is then initialized by
  fCdf->InitialiseForComponent(Z, numEnergies);
  
  // Now the data files are read in for all energies. At each energy, 
  // there are fNumAngles angles and fNumAngles CDF values.

  std::vector<G4PhysicsVector*> tempDataCDF;
  // Loop through each energy
  for (int j=0; j<numEnergies; j++)
    {
      // Push back a new G4PhysicsVector for the jth energy
      tempDataCDF.push_back(new G4PhysicsVector(false));
          
      // For use with David's PhysicsVector class
      //tempDataCDF.push_back(new G4PhysicsVector(false,false,true));
          
      //tempDataCDF.push_back(new G4PhysicsVector(false,false,false));
      //tempDataCDF.push_back(new G4PhysicsVector());
    } 
          
  // Open a temporary stream. The data for the jth energy group is copied 
  // to the temp file and then the temp file is sent to 
  // G4PhysicsVector::Retrieve(). Once the data is stored in the 
  // G4PhysicsVector, tempDataCDF, it is then passed to the ElementData, cdf, 
  // which is the container for all of the data for each element and energy.
  
  //static G4Mutex m = G4MUTEX_INITIALIZER;
  //G4AutoLock l(&m);
  
  // Open the stream and call file "tempDataFile.dat"
  std::ofstream file("tempDataFile.dat",
                     std::fstream::out | std::fstream::trunc);
  
  // Write to file the lower/upper bounds and the number of grid points 
  // for each column of data. The data in tempfile is 2xfNumAngles, 
  // hence fNumAngles fNumAngles. 
  // file<<lowerBound<<" "<<upperBound<<" "<<fNumAngles<<" "
  // <<fNumAngles<<G4endl;
  file<<"-1. 1. "<<fNumAngles<<" "<<fNumAngles<<G4endl;

  // Start while loop over the entire data file opened above 
  // e.g. pwe_cdf_79.dat. This data contains all of the data for each energy 
  G4int cntr = 0;
  G4int j = 0;
  G4double temp1, temp2;
  while(in2>>temp1>>temp2)
  {        
          // Write data to temporary file
    file<<setprecision(16)<<temp1<<"  "<<temp2<<G4endl;  
    cntr++;
    
    // When the first fNumAngles data points are copied to tempDataFile.dat, 
    // store data in tempDataCDF[j]. Then increment the energy index, j, 
    // close and clear the streams, and then reopen the streams such that 
    // the next fNumAngles data points are copied to tempDataFile.dat.
    if (cntr==fNumAngles) 
    {
      cntr=0;
      std::ifstream inTemp("tempDataFile.dat", 
                           std::ifstream::binary|std::ifstream::in);
      tempDataCDF[j]->Retrieve(inTemp,true);
      fCdf->AddComponent(Z,j,tempDataCDF[j]); 
      j++;
      inTemp.close(), inTemp.clear(), file.close(), file.clear();
      file.open("tempDataFile.dat",std::fstream::out | std::fstream::trunc);
      file<<"-1. 1. "<<fNumAngles<<" "<<fNumAngles<<G4endl;
    }
  }
    
  in2.close();
  
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DiscreteScatteringModel::SampleSecondaries(
                               std::vector<G4DynamicParticle*>*,
                               const G4MaterialCutsCouple* couple,
                               const G4DynamicParticle* p,
                               G4double cutEnergy, 
                               G4double)
{
  G4double E = p->GetKineticEnergy();

  if(E < fLowEnergyLimit) return;
        
  // Select random atom
  G4int Z = G4lrint(SelectRandomAtom(couple,p->GetDefinition(),
                                     E,cutEnergy,E)->GetZ());

  // Determine the energy bin
  G4double logE = G4Log(E);
  G4ThreeVector dir = p->GetMomentumDirection(); //old direction
  G4int         j   = fTcs->GetElementData(Z)->FindBin(logE,0); 

  //-------------------------------------------------------------------------
  // it would be nice to have the following block of code in G4PhysicsVector.
  // it could be a simple function. 
  
  //    v.GetVectorLength()  - number of points
  //    v.Energy(size_t idx) - value x(i)
  //    v[i]                 - value y(i)
  
  // This is a monte carlo interpolation scheme
  G4double e1 = fTcs->GetElementData(Z)->Energy(j);
  G4double e2 = fTcs->GetElementData(Z)->Energy(j+1);
  
  // This is a monte carlo interpolation scheme
  G4double pie1 = (logE-e1)/(e2-e1);
  G4double    r = G4UniformRand();
  if (r<pie1){ ++j; } 
  
  //-------------------------------------------------------------------------
  // it would be nice to have the following block of code in G4PhysicsVector.
  //---------------
  // Given the energy grid value associated with the DCS, 
  // sample a deflection cosine       
  r       = G4UniformRand();
  G4int k = -1;
  // First test if the angle is the most probable angle
  if (r>(*fCdf->GetComponentDataByIndex(Z,j)).Energy(fNumAngles-2)) 
    { k = fNumAngles-1; }
  // or the least probable... not sure why I do this (maybe because 
  // it is a simple check)
  else if (r<=(*fCdf->GetComponentDataByIndex(Z,j)).Energy(0)) { k = 0; }
  // if neither then loop through remaining angles, break when locating angle
  else {
    for (G4int i=fNumAngles-2; i>0; --i) { 
      if ( (r>(*fCdf->GetComponentDataByIndex(Z,j)).Energy(i-1)) 
           && (r<= (*fCdf->GetComponentDataByIndex(Z,j)).Energy(i)) )
        { k=i; break;}
    }
  }
  // Throw an error if an angle was not sampled, data is probably no good  
  if(k<0)
  {
    G4cout << "G4DiscreteScatteringModel::SampleSecondaries():"
           << " CDF was not inverted properly "<<k<<G4endl;
    for (G4int i=0;i<fNumAngles;++i) { 
      G4cout<<i<<" "<<(*fCdf->GetComponentDataByIndex(Z,j)).Energy(i)<<" "<<r<<G4endl;
    }
  }
  //-------------------------------------------------------------------------

  // Otherwise, go get the angle and pass it too local method GetNewDirection.
  // Then do transformation and update fParticleChange.
  
  //G4cout<<"------------------"<<G4endl;
  //G4cout<<"Sample Secondaries"<<G4endl;
  //G4cout<<G4Exp(e1)<<" "<<E<<" "<<G4Exp(e2)<<" "<<pie1<<" "
  // <<(*fCdf->GetComponentDataByIndex(Z,j)).Energy(k)<<G4endl;
  //G4cout<<"------------------"<<G4endl;
  //G4cout<<" "<<G4endl;    
  
  G4ThreeVector newDirection = 
    GetNewDirection((*fCdf->GetComponentDataByIndex(Z,j))[k]);
  newDirection.rotateUz(dir);   
  fParticleChange->ProposeMomentumDirection(newDirection);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
                
G4double G4DiscreteScatteringModel::ComputeCrossSectionPerAtom(
  const G4ParticleDefinition*,G4double E, G4double Z, G4double,
                              G4double, G4double)
{
  // super simple, very compact look-up for total cross section 
  if (E==0.) {return 1e30;}
  //G4cout<<"------------------"<<G4endl;
  //G4cout<<"ComputeCrossSectionPerAtom"<<G4endl;
  //G4cout<< E<<" "<<G4Log(E)<<" "<< G4Exp(fTcs->GetValueForElement(Z,log(E)))
  //<<" "<<fTcs->GetValueForElement(Z,log(E))<<G4endl;
  //G4cout<<"------------------"<<G4endl;
  //G4cout<<" "<<G4endl;  
  
  //G4cout<< E<<" "<<G4Exp(fTcs->GetValueForElement(Z,log(E)))<<G4endl;
  
  return G4Exp(fTcs->GetValueForElement(Z,G4Log(E))); 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector G4DiscreteScatteringModel::GetNewDirection(G4double z1)
{
  G4ThreeVector dir(0.0,0.0,1.0);

  G4double sint = sin(acos(z1));
  G4double cost = sqrt(1.0 - sint*sint);
  G4double phi  = twopi* G4UniformRand();
  G4double dirx = sint*cos(phi);
  G4double diry = sint*sin(phi);
  G4double dirz = cost;

  dir.set(dirx,diry,dirz);
  return dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


