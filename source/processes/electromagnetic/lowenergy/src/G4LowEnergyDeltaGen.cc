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
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      ------------ G4LowEnergyDeltaGen --------------  
//
// **************************************************************
// 
// Modified: 03.09.01 First implementation.
//
// --------------------------------------------------------------

#include "G4LowEnergyDeltaGen.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4LogLogInterpolation.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4DataVector.hh"
#include "G4ShellEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh" 
#include "g4std/fstream"
#include "g4std/strstream"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4LowEnergyDeltaGen::G4LowEnergyDeltaGen(
  const G4String& nam) : G4VEMSecondaryGenerator(nam),   
  length_a(16),
  length_p(14),
  length_z(99),
  interpolation(0),
  interpolation1(0),
  bindingData(0),
  lowEnergyLimit(250.0*eV),
  highEnergyLimit(10.0*GeV)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4LowEnergyDeltaGen::~G4LowEnergyDeltaGen()
{
  Clear();
  if(interpolation)  delete interpolation;
  if(interpolation1) delete interpolation1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyDeltaGen::Clear()
{
  // Reset the map of data sets: remove the data sets from the map 
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::iterator pos;

  for (pos = param.begin(); pos != param.end(); pos++)
    {
      G4VEMDataSet* dataSet = pos->second;
      param.erase(pos);
      delete dataSet;
    }

  G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::iterator vec;

  for (vec = energyLimit1.begin(); vec != energyLimit1.end(); vec++)
    {
      G4DataVector* data = vec->second;
      data->clear();
      energyLimit1.erase(vec);
      delete data;
    }

  for (vec = energyLimit2.begin(); vec != energyLimit2.end(); vec++)
    {
      G4DataVector* data = vec->second;
      data->clear();
      energyLimit2.erase(vec);
      delete data;
    }

  activeZ.clear();
  if(bindingData) delete bindingData;
  if(verbose > 0) {
    G4cout << "G4LowEnergyDeltaGen is cleared"
           << G4endl;
      }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyDeltaGen::Initialize()
{
  if(verbose > 0) {
    G4cout << "G4LowEnergyDeltaGen::Initialize start"
           << G4endl;
      }

  if(!interpolation1) interpolation1 = new G4LogLogInterpolation();
  if(!interpolation)  interpolation  = new G4SemiLogInterpolation();

  if(verbose > 1) G4cout << "G4LogLogInterpolation is used" << G4endl;
  if(verbose > 1) G4cout << "G4SemiLogInterpolation is used" << G4endl;

  // define active elements

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
     G4Exception("G4CrossSectionHandler: no MaterialTable found)");

  G4int nMaterials = materialTable->length();
  
  for (G4int m=0; m<nMaterials; m++) {

    const G4Material* material= (*materialTable)[m];        
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4int nElements = material->GetNumberOfElements();
      
    for (size_t iEl=0; iEl<nElements; iEl++) {
      G4Element* element = (*elementVector)[iEl];
      G4double Z = element->GetZ();
      if (!(activeZ.contains(Z)) && Z > 0 && Z < length_z) {
	    activeZ.push_back(Z);
            if(verbose > 1) G4cout << "New active Z= " << Z << G4endl;
      }
    }
  }

  // Read parameters 
  bindingData = new G4ShellEMDataSet(0, "fluor/binding.dat", interpolation1, 
                                     1., 1.);
  if(verbose > 1) G4cout << "Binding data have been read" << G4endl;

   
  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep = "G4EMDataSet - G4LEDATA environment variable not set";
      G4Exception(excep);
    }
  G4String pathString(path);
  pathString += "/ioni/io-co-";  

  G4double energy;
  G4DataVector* e;
  G4VEMDataSet* set;
  G4std::vector<G4DataVector*> a;
  G4std::vector<G4ShellEMDataSet*> p;
  a.resize(length_a);
  p.resize(length_p);

  size_t nZ = activeZ.size();
  G4DataVector* tt1;
  G4DataVector* tt2;

  if(verbose > 1) G4cout << "Start read shell parameters nZ= " << nZ << G4endl;

  for (G4int i=0; i<nZ; i++) {

    tt1 = new G4DataVector();
    tt2 = new G4DataVector();
    
    G4int Z = (G4int) activeZ[i];      
    char nameChar[100] = {""};
    G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
    ost << pathString << Z << ".dat";
    G4String name(nameChar);  

    G4std::ifstream file(name);
    G4std::filebuf* lsdp = file.rdbuf();
  
    if (! (lsdp->is_open()) ) {
      G4String excep = G4String("G4LowEnergyDeltaGen::Initialize: ")
                     + name;
      G4Exception(excep);
    }
    if(verbose > 1) {
      G4cout << "G4LowEnergyDeltaGen: the file <"
             << name << "> is opened for reading"
             << G4endl;
    }

    // The file is organized into two columns:
    // 1st column is the energy
    // 2nd column is the corresponding value
    // The file terminates with the pattern: -1   -1
    //                                       -2   

    e   = new G4DataVector();
    for(size_t j=0; j<length_a; j++) {
      a[j] = new G4DataVector();
    } 
    for(size_t k=0; k<length_p; k++) {
      G4int id = i*20 + k;
      p[k] = new G4ShellEMDataSet(id, interpolation);
    }

    if(verbose > 2) G4cout << "Start read the file i= " << i << G4endl;

    do {
      file >> energy;
      if (energy == -2) break;
      if (energy >  -1) e->push_back(energy);

      G4double x;
      for (size_t j=0; j<length_a; j++) {
          file >> x;    
          if(energy > -1) a[j]->push_back(x);
      }

      // fill map
      if(energy < 0) {

        G4double tcut1 = (*a[14])[0];
        G4double tcut2 = (*a[15])[0];
        G4int id;
        for(size_t k=0; k<14; k++) {
            id = Z*20 + k;
            set = new G4EMDataSet(id, e, a[k], interpolation);
            p[k]->AddComponent(set);
        } 
        if(tcut1 == 0.0) tcut1 = DBL_MAX;
        if(tcut2 == 0.0) tcut2 = DBL_MAX;
        tt1->push_back(tcut1);
        tt2->push_back(tcut2);
      }
    } while (energy > -2);

    file.close();

    energyLimit1[Z] = tt1;
    energyLimit2[Z] = tt2;

    for(G4int kk=0; kk<length_p; kk++) {
      G4int id = Z*20 + kk;
      param[id] = p[kk];
      if(verbose > 2) {
        G4cout << "G4LowEnergyDeltaGen: parameters ID= " << id
               << " are saved" 
               << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyDeltaGen::Probability(G4int atomicNumber,
                                          G4int shellNumber,
	  	                          G4double kineticEnergy,
                                          G4double tcut,
                                          G4double thigh) const
{
  // Control limits
  G4double tmax = thigh;
  if(tmax > kineticEnergy) tmax = MaxSecondaryEnergy(G4Electron::Electron(), 
                                                     kineticEnergy);
  if(tcut >= tmax) return 0.0;

  // Access parameters
  const G4DataVector* p = FindParameters(atomicNumber, shellNumber, 
                                        kineticEnergy);
  G4double t1 = FindEnergyLimit1(atomicNumber, shellNumber);
  G4double t2 = FindEnergyLimit2(atomicNumber, shellNumber);

  G4double z = (G4double) atomicNumber;
  const G4VEMDataSet* data = bindingData->GetComponent(shellNumber);
  G4double bindingEnergy = data->FindValue(z);
    
  G4double val = 0.0;
  G4double nor = 0.0;
  
  // Integration is performed in 3 energy arias, in each probability 
  // is integrated between a2 and a3, normalisation - between a1 and a3

  // --- First function ---

  G4double a1 = 0.0; 
  G4double a2 = G4std::min(tcut,t1); 
  G4double a3 = G4std::min(tmax,t1);

  if(a2 < a3) val += IntSpectrum1(0, 0, a2, a3, bindingEnergy, (*p));
  nor += IntSpectrum1(0, 0, a1, a3, bindingEnergy, (*p));

  // --- Second function ---

  if(tmax > t1) {
    G4double a1 = t1; 
    G4double a2 = G4std::min(G4std::max(tcut,t1),t2); 
    G4double a3 = G4std::min(t1,t2);

    if(a2 < a3) val += IntSpectrum2(0, a2, a3, (*p));
    if(a1 < a3) nor += IntSpectrum2(0, a1, a3, (*p));
  }

  // --- Third function ---

  if(tmax > t2) {

    G4double a1 = t2; 
    G4double a2 = G4std::max(tcut,t2);
    G4double a3 = tmax;

    if(a2 < a3) val += IntSpectrum1(0, 1, a2, a3, bindingEnergy, (*p));
    if(a1 < a3) nor += IntSpectrum1(0, 1, a1, a3, bindingEnergy, (*p));
  }

  if(nor > 0.0) val /= nor;
  else          val = 0.0;
  delete p;

  return val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyDeltaGen::AverageEnergy(G4int atomicNumber,
                                            G4int shellNumber,
	  	                            G4double kineticEnergy,
                                            G4double tcut) const
{
  // Control limits
  G4double tmax = G4std::min(tcut, MaxSecondaryEnergy(G4Electron::Electron(), 
                                                      kineticEnergy));
  if(verbose > 1) {
    G4cout << "G4LowEnergyDeltaGen::AverageEnergy: Z= " << atomicNumber
           << "; shell= " << shellNumber
           << "; E(keV)= " << kineticEnergy/keV
           << G4endl;
  }

  // Access parameters
  const G4DataVector* p = FindParameters(atomicNumber, shellNumber, 
                                         kineticEnergy);
  G4double t1 = FindEnergyLimit1(atomicNumber, shellNumber);
  if(atomicNumber == 14) G4cout << "t1= " << t1 << G4endl;
  G4double t2 = FindEnergyLimit2(atomicNumber, shellNumber);
  if(atomicNumber == 14) G4cout << "t2= " << t2 << G4endl;
    
  G4double z = (G4double) atomicNumber;
  const G4VEMDataSet* data = bindingData->GetComponent(shellNumber);
  if(atomicNumber == 14) data->PrintData();
  G4double bindingEnergy = data->FindValue(z);
  if(verbose > 1) {
    G4cout << "G4LowEnergyDeltaGen::AverageEnergy: Z= " << atomicNumber
           << "; shell= " << shellNumber
           << "; bindingE(keV)= " << bindingEnergy/keV
           << G4endl;
  }

  G4double val = 0.0;
  G4double nor = 0.0;
  
  // Integration is performed in 3 energy arias, in each energy and 
  // normalisation is integrated between a1 and a2

  // --- First function ---

  G4double a1 = 0.0; 
  G4double a2 = G4std::min(tmax,t1);

  val += IntSpectrum1(1, 0, a1, a2, bindingEnergy, (*p));
  nor += IntSpectrum1(0, 0, a1, a2, bindingEnergy, (*p));
  if(verbose > 1) {
    G4cout << "1-st function val= " << val << "; nor= " << nor << G4endl;
  }

  // --- Second function ---

  if(tmax > t1) {
    G4double a1 = t1; 
    G4double a2 = G4std::min(tmax,t2); 

    if(a1 < a2) {
      val += IntSpectrum2(1, a1, a2, (*p));
      nor += IntSpectrum2(0, a1, a2, (*p));
    }
  }
  if(verbose > 1) {
    G4cout << "2-nd function val= " << val << "; nor= " << nor << G4endl;
  }

  // --- Third function ---

  if(tmax > t2) {

    G4double a1 = t2; 
    G4double a2 = G4std::max(tmax,t2);

    if(a1 < a2) {
      val += IntSpectrum1(1, 1, a1, a2, bindingEnergy, (*p));
      nor += IntSpectrum1(0, 1, a1, a2, bindingEnergy, (*p));
    }
  }
  if(verbose > 1) {
    G4cout << "3-d  function val= " << val << "; nor= " << nor << G4endl;
  }

  if(nor > 0.0) val /= nor;
  else          val = 0.5*tcut;
  delete p;

  return val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyDeltaGen::GenerateSecondary(
                             const G4DynamicParticle* aParticle,
			           G4ParticleChange* theChange,
                                   G4int atomicNumber, 
                                   G4int shellNumber,
                                   G4double tcut,
                                   G4double thigh)
{
  // Delta electron production mechanism on base of the model
  // J. Stepanek " A program to determine the radiation spectra due 
  // to a single atomic subshell ionisation by a particle or due to 
  // deexcitation or decay of radionuclides", 
  // Comp. Phys. Comm. 1206 pp 1-19 (1997)

  G4double kineticEnergy = aParticle->GetKineticEnergy();
  const G4ParticleDefinition* particle = aParticle->GetDefinition();

  // Control limits
  G4double tmax = thigh;
  if(tmax > kineticEnergy) tmax = MaxSecondaryEnergy(particle, kineticEnergy);
  if(tcut >= tmax) return;

  // Access parameters
  const G4DataVector* p = FindParameters(atomicNumber, shellNumber, 
                                         kineticEnergy);
  G4double t1 = FindEnergyLimit1(atomicNumber, shellNumber);
  G4double t2 = FindEnergyLimit2(atomicNumber, shellNumber);
    
  G4double z = (G4double) atomicNumber;
  const G4VEMDataSet* data = bindingData->GetComponent(shellNumber);
  G4double bindingEnergy = data->FindValue(z);

  G4double fmax = G4std::max((*p)[6], (*p)[11]);

  // Sample delta energy

  G4double tdel;
  G4double f, q;
  do {
    tdel = tcut + G4UniformRand()*(tmax - tcut);
    if(tdel < t1)      f = Spectrum1(0, tdel, bindingEnergy, (*p));
    else if(tdel < t2) f = Spectrum2(tdel, (*p));
    else               f = Spectrum1(1, tdel, bindingEnergy, (*p));
    q = fmax * G4UniformRand();
    if(f > fmax) {
      G4cout << "G4LowEnergyDeltaGen::GenerateSecondary WARNING: "
             << " f= " << f << " > " << " fmax= " << fmax
             << " Z= " << atomicNumber
             << " shell= " << shellNumber
             << " E(MeV)= " << kineticEnergy/MeV
             << G4endl;
    }
  } while (q > f);



  // Transform to shell potential
  G4double deltaKinE = tdel + 2.0*bindingEnergy; 
  G4double primaryKinE = kineticEnergy + 2.0*bindingEnergy;   

  // sampling of scattering angle neglecting atomic motion
  G4double deltaMom = sqrt(deltaKinE*(deltaKinE + 2.0*electron_mass_c2));
  G4double primaryMom = sqrt(primaryKinE*(primaryKinE + 2.0*electron_mass_c2));
   
  G4double cost = deltaKinE * (primaryKinE + 2.0*electron_mass_c2)
                            / (deltaMom * primaryMom);

  if (cost > 1.) cost = 1.;
  G4double sint = sqrt(1. - cost*cost);
  G4double phi  = twopi * G4UniformRand(); 
  G4double dirx = sint * cos(phi);
  G4double diry = sint * sin(phi);
  G4double dirz = cost;

  // Rotate to incident electron direction
  G4ThreeVector primaryDirection = aParticle->GetMomentumDirection();
  G4ThreeVector deltaDir(dirx,diry,dirz);
  deltaDir.rotateUz(primaryDirection);
  dirx = deltaDir.x();
  diry = deltaDir.y();
  dirz = deltaDir.z();

  // Take into account atomic motion del is relative momentum of the motion
  // kinetic energy of the motion == bindingEnergy in V.Ivanchenko model

  cost = 2.0*G4UniformRand() - 1.0;
  sint = sqrt(1. - cost*cost);
  phi  = twopi * G4UniformRand(); 
  G4double del = sqrt(bindingEnergy *(bindingEnergy + 2.0*electron_mass_c2))
               / deltaMom;
  dirx += del* sint * cos(phi);
  diry += del* sint * sin(phi);
  dirz += del* cost;
  G4double norm = 1.0/sqrt(dirx*dirx + diry*diry + dirz*dirz);
  dirx *= norm;
  diry *= norm;
  dirz *= norm;

  // Find out new primary electron direction
  G4double finalPx = primaryMom*primaryDirection.x() - deltaMom*dirx; 
  G4double finalPy = primaryMom*primaryDirection.y() - deltaMom*diry; 
  G4double finalPz = primaryMom*primaryDirection.z() - deltaMom*dirz; 

  // create G4DynamicParticle object for delta ray
  theChange->SetNumberOfSecondaries(1);
  G4DynamicParticle* theDeltaRay = new G4DynamicParticle();
  theDeltaRay->SetKineticEnergy(tdel);
  theDeltaRay->SetMomentumDirection(dirx, diry, dirz); 
  theDeltaRay->SetDefinition(G4Electron::Electron());
  theChange->AddSecondary(theDeltaRay);
     
  G4double theEnergyDeposit = bindingEnergy;

  // Fluorescence should be implemented here


  // fill ParticleChange 
  // changed energy and momentum of the actual particle

  G4double finalKinEnergy = kineticEnergy - tdel - theEnergyDeposit;
  if(finalKinEnergy < 0.0) {
    theEnergyDeposit += finalKinEnergy;
    finalKinEnergy    = 0.0;
  }
  norm = 1.0/sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz);
  finalPx *= norm;
  finalPy *= norm;
  finalPz *= norm;

  theChange->SetMomentumChange(finalPx, finalPy, finalPz);
  theChange->SetEnergyChange(finalKinEnergy);
  if(theEnergyDeposit < 0.) {
    G4cout << "G4LowEnergyIonisation: Negative energy deposit: " 
           << theEnergyDeposit/eV << " eV" << G4endl;
    theEnergyDeposit = 0.0;
  }
  theChange->SetLocalEnergyDeposit(theEnergyDeposit);
  delete p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4DataVector* G4LowEnergyDeltaGen::FindParameters(G4int atomicNumber, 
                                                        G4int shellNumber, 
                                                        G4double e) const
{
  G4DataVector* par = new G4DataVector();
  G4int id = atomicNumber*20;
  
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;

  for(G4int i=0; i<length_p; i++) { 

    G4int index = id + i;
    pos = param.find(index);
    if (pos!= param.end()) {
      G4VEMDataSet* dataSet = pos->second;
      size_t nShells = dataSet->NumberOfComponents();

      if(shellNumber < nShells) { 
        const G4VEMDataSet* comp = dataSet->GetComponent(shellNumber);
        const G4DataVector ener = comp->GetEnergies(0);
        G4double ee = G4std::max(ener.front(),G4std::min(ener.back(),e));
        G4double value = comp->FindValue(ee);
        if(verbose > 2) {
          G4cout << "FindParameters: id= " << id 
                 << "; e= " << e
                 << "; ee= " << ee
                 << "; p= " << value
                 << G4endl;
	}
        par->push_back(value);

      } else {
        G4cout << "WARNING: G4LowEnergyDeltaGen::FindParameters "
               << "has no parameters for shell " << shellNumber 
               << " for Z= " << atomicNumber 
               << G4endl;
      }
    } else {
      G4cout << "WARNING: G4LowEnergyDeltaGen::FindValue "
             << "did not find ID = "
	     << index << G4endl;
    }
  }
  return par;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyDeltaGen::FindEnergyLimit1(G4int atomicNumber, 
                                               G4int shellNumber) const
{
  G4double value = 0.0;
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::const_iterator pos;
  pos = energyLimit1.find(atomicNumber);

  if (pos!= param.end()) {
      G4DataVector* data = pos->second;
      G4int nShells = data->size();
      if(shellNumber < nShells) 
	{ 
          value = (*data)[shellNumber];
        }
      else
        {
          G4cout << "WARNING: G4LowEnergyDeltaGen::FindEnergyLimit1 "
                 << "has no parameters for shell "
	         << shellNumber 
                 << " for Z= "
                 << atomicNumber << G4endl;
	}
    }
  else
    {
      G4cout << "WARNING: G4LowEnergyDeltaGen::FindEnergyLimit1 "
             << "did not find Z = "
	     << atomicNumber << G4endl;
    }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyDeltaGen::FindEnergyLimit2(G4int atomicNumber, 
                                               G4int shellNumber) const
{
  G4double value = 0.0;
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::const_iterator pos;
  pos = energyLimit2.find(atomicNumber);

  if (pos!= param.end()) {
      G4DataVector* data = pos->second;
      G4int nShells = data->size();
      if(shellNumber < nShells) 
	{ 
          value = (*data)[shellNumber];
        }
      else
        {
          G4cout << "WARNING: G4LowEnergyDeltaGen::FindEnergyLimit2 "
                 << "has no parameters for shell "
	         << shellNumber 
                 << " for Z= "
                 << atomicNumber << G4endl;
	}
    }
  else
    {
      G4cout << "WARNING: G4LowEnergyDeltaGen::FindEnergyLimit2 "
             << "did not find Z = "
	     << atomicNumber << G4endl;
    }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyDeltaGen::IntSpectrum1(size_t moment, 
                                           size_t function, 
                                           G4double tmin,
                                           G4double tmax,
                                           G4double bindingEnergy,
                                     const G4DataVector& p) const
{
  G4double value = 0.0;
  if (tmin >= tmax) return value;
  size_t i = 0;
  size_t imax = 6;
  if(function) {
    i = 7;
    imax = 4;
  }
  for(size_t j=0; j<imax; j++) {
    value += IntSpectrum3(moment, j+2, tmin, tmax, bindingEnergy) * p[j + i];
  } 
  if(verbose > 1) {
    G4cout << "moment= " << moment << "; i= " << i 
           << "; imax= " << imax
           << "; tmin= " << tmin
           << "; tmax= " << tmax
           << "; b= " << bindingEnergy
           << "; p= " << p[i]
           << "; Int= " << value
           << G4endl;
  }

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyDeltaGen::IntSpectrum2(size_t moment, 
                                           G4double tmin,
                                           G4double tmax,
                                     const G4DataVector& p) const
{
  G4double value = 0.0;
  if (tmin >= tmax) return value;
  G4double c2 = p[13];

  // Integral on 1-st momentum
  if(moment) {
    value = (pow(tmax, c2 + 2.) - pow(tmin, c2 + 2.)) / (c2 + 2.); 

  // Integral on 0-st momentum
  } else {
    value = (pow(tmax, c2 + 1.) - pow(tmin, c2 + 1.)) / (c2 + 1.); 
  } 
  value *= p[12];
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyDeltaGen::IntSpectrum3(size_t moment, 
                                           size_t power, 
                                           G4double tmin,
                                           G4double tmax,
                                           G4double b) const
{
  G4double value = 0.0;
  if (tmin >= tmax) return value;
  G4double y = 1.0 - (G4double)power;
 
  value = (pow(tmax + b, y) - pow(tmin + b, y)) / y; 

  // Integral on 1-st momentum
  if(moment) {
    
    value *= (-b);
    y += 1.0;

    if(y < -1.0) value += (pow(tmax + b, y) - pow(tmin + b, y)) / y;

    else         value += log((tmax + b)/(tmin + b));

  } 

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyDeltaGen::Spectrum1(size_t function, 
                                        G4double tdelta,
                                        G4double bindingEnergy,
                                  const G4DataVector& p) const
{
  G4double value = 0.0;
  size_t i = 0;
  size_t imax = 6;
  if(function) {
    i = 7;
    imax = 4;
  }
  for(size_t j=0; j<imax; j++) {
    G4double y = -2.0 - (G4double)j;
    value += pow(tdelta + bindingEnergy, y) * p[j + i];
  } 

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyDeltaGen::Spectrum2(G4double tdelta,
                                  const G4DataVector& p) const
{
  G4double value = p[12] * pow(tdelta, p[13]);
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyDeltaGen::PrintData() const
{
  G4cout << generatorName 
         << ": e- Delta generation using EEDL database, "
         << "\n       gamma energy sampled from a parametrised formula."
         << "\n       Energy range = ["
         << lowEnergyLimit/eV << "eV, "
         << highEnergyLimit/GeV << "GeV]" 
         << G4endl;
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....












