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
//      ------------ G4LowEnergyBremsstrahlungGen --------------  
//
// **************************************************************
// 
// Modified: 03.09.01 First implementation.
//
// --------------------------------------------------------------

#include "G4LowEnergyBremsstrahlungGen.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh" 
#include "g4std/fstream"
#include "g4std/strstream"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4LowEnergyBremsstrahlungGen::G4LowEnergyBremsstrahlungGen(
  const G4String& nam) : G4VEMSecondaryGenerator(nam),   
  length(99),
  interpolation(0),
  lowestEnergyGamma(0.1*eV),
  lowEnergyLimit(250.0*eV),
  highEnergyLimit(10.0*GeV)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4LowEnergyBremsstrahlungGen::~G4LowEnergyBremsstrahlungGen()
{
  if(interpolation) delete interpolation;
  Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyBremsstrahlungGen::Clear()
{
  // Reset the map of data sets: remove the data sets from the map 
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::iterator pos;

  for (pos = paramA.begin(); pos != paramA.end(); pos++)
    {
      G4VEMDataSet* dataSet = pos->second;
      paramA.erase(pos);
      delete dataSet;
    }

  c.clear();
  d.clear();
  activeZ.clear();
  delete interpolation;
  if(verbose > 0) {
    G4cout << "G4LowEnergyBremsstrahlungGen is cleared"
           << G4endl;
      }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyBremsstrahlungGen::Initialize()
{
  if(verbose > 0) {
    G4cout << "G4LowEnergyBremsstrahlungGen::Initialize start"
           << G4endl;
      }

  if(!interpolation) interpolation = new G4LogLogInterpolation();

  if(verbose > 1) G4cout << "G4LogLogInterpolation is used" << G4endl;

  // define active elements

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
     G4Exception("G4CrossSectionHandler: no MaterialTable found)");

  G4int nMaterials = materialTable->length();
  
  for (G4int m=0; m<nMaterials; m++) {

    const G4Material* material= (*materialTable)[m];        
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4int nElements = material->GetNumberOfElements();
      
    for (G4int iEl=0; iEl<nElements; iEl++) {
      G4Element* element = (*elementVector)[iEl];
      G4double Z = element->GetZ();
      if (!(activeZ.contains(Z)) && Z > 0 && Z < length) {
	    activeZ.push_back(Z);
            if(verbose > 1) G4cout << "New active Z= " << Z << G4endl;
      }
    }
  }


  // Read parameters A
  if(verbose > 0) {
    G4cout << "G4LowEnergyBremsstrahlungGen read G4LEDATA/brem/br-co-a.dat"
           << G4endl;
      }
   
  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep = "G4EMDataSet - G4LEDATA environment variable not set";
      G4Exception(excep);
    }

  G4String pathString(path);
  G4String name_a = pathString + "/brem/br-co-a.dat";  
  G4std::ifstream file_a(name_a);
  G4std::filebuf* lsdp_a = file_a.rdbuf();
  
  if (! (lsdp_a->is_open()) ) {
     G4String excep = G4String("G4LowEnergyBremsstrahlungGen::Initialize: ")
                    + name_a;
     G4Exception(excep);
  }

  // The file is organized into two columns:
  // 1st column is the energy
  // 2nd column is the corresponding value
  // The file terminates with the pattern: -1   -1
  //                                       -2   -2

  G4DataVector* energies;
  G4DataVector* data;
  G4double e = 0;
  G4double a = 0;
  size_t nZ  = activeZ.size();
  energies   = new G4DataVector(); 
  data       = new G4DataVector(); 
  G4int z    = 0;

  do {
    file_a >> e >> a;

    // End of the next element
    if (e == -1) {

      z++;
      G4bool used = false;
      for (size_t i=0; i<nZ; i++) {
    
	// fill map
        if(z == (G4int) activeZ[i]) {
          G4VEMDataSet* dataSet = 
                        new G4EMDataSet(z, energies, data, interpolation);
          paramA[z] = dataSet;
          used = true;
	}
      }
      if(!used) {
        energies->clear();
        data->clear();
      }
  
      if (z > 98) break;

      energies   = new G4DataVector(); 
      data       = new G4DataVector(); 
      
    } else {

      energies->push_back(e);
      data->push_back(a);
    }
  } while (e != -2);
  
  file_a.close();

  // Read parameters B
  if(verbose > 0) {
    G4cout << "G4LowEnergyBremsstrahlungGen read G4LEDATA/brem/br-co-b.dat"
           << G4endl;
      }

  G4String name_b = pathString + "/brem/br-co-b.dat";  
  G4std::ifstream file_b(name_b);
  G4std::filebuf* lsdp_b = file_b.rdbuf();
  
  if (! (lsdp_b->is_open()) ) {
     G4String excep = G4String("G4LowEnergyBremsstrahlungGen::InitializeMe: ")
                    + name_b;
     G4Exception(excep);
  }

  // The file is organized into two columns:
  // 1st column is c; 2nd column is d;
  c.clear();
  d.clear();

  for (size_t j=0; j<length; j++) {

    file_b >> e >> a;
    if(e == -1) break;
    c.push_back(e);
    d.push_back(a);
  }
  file_b.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyBremsstrahlungGen::Probability(G4int atomicNumber,
                                                   G4int,
	  	                                   G4double kineticEnergy,
                                                   G4double tmin,
                                                   G4double tmax) const
{
  // Control limits
  G4double emin = G4std::max(tmin, lowestEnergyGamma);
  G4double emax = G4std::min(tmax, kineticEnergy);
  
  if(emin >= emax) return 0.0;

  G4double a = FindValueA(atomicNumber, kineticEnergy);
  size_t   i = atomicNumber - 1;
  G4double b = c[i]*log10(kineticEnergy/MeV) + d[i];
  if(a < 0) a = 0.0;
  if(b < 0) b = 0.0;
    
  G4double val = a*log(emax/emin) + b*(emax - emin);
  G4double nor = a*log(emax/lowestEnergyGamma) + b*(emax - lowestEnergyGamma);
  return val/nor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyBremsstrahlungGen::AverageEnergy(G4int atomicNumber,
                                                     G4int,
	  	                                     G4double kineticEnergy,
                                                     G4double tcut) const
{
  // Control limits
  G4double emin = lowestEnergyGamma;
  G4double emax = G4std::min(tcut, kineticEnergy);
  
  if(emin >= emax) return 0.0;

  G4double a = FindValueA(atomicNumber, kineticEnergy);
  size_t   i = atomicNumber - 1;
  G4double b = c[i]*log10(kineticEnergy/MeV) + d[i];
  if(a < 0) a = 0.0;
  if(b < 0) b = 0.0;
    
  G4double val = a*(emax-emin) + 0.5*b*(emax*emax - emin*emin);
  G4double nor = a*log(emax/emin) + b*(emax - emin);
  return val/nor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyBremsstrahlungGen::GenerateSecondary(
                             const G4DynamicParticle* aParticle,
			           G4ParticleChange* theChange,
                                   G4int atomicNumber, 
                                   G4int,
                                   G4double tmin,
                                   G4double tmax)
{
  // This parametrization is derived from : 
  // Migdal corrections (dielectric suppression). 
  // Migdal: Phys Rev 103:1811 (1956); 
  // Messel & Crawford: Pergamon Press (1970)
  
  theChange->SetNumberOfSecondaries(1);
  G4double kineticEnergy = aParticle->GetKineticEnergy();

  // Control limits
  G4double emin = G4std::max(tmin, lowestEnergyGamma);
  G4double emax = G4std::min(tmax, kineticEnergy);
  if(emin >= emax) return;

  G4double a = FindValueA(atomicNumber, kineticEnergy);
  size_t   i = atomicNumber - 1;
  G4double b = c[i]*log10(kineticEnergy/MeV) + d[i];
  if(a < 0) a = 0.0;
  if(b < 0) b = 0.0;

  // Sample gamma energy
  // The emitted gamma energy is from EEDL data fitted with A/E+B function.
  // Original formula A/E+B+C*E and sampling methods are reported by 
  // J.Stepanek formula has been modified by A. Forti and S. Giani. 

  G4double amaj = a + b*emax;
  G4double amax = log(emax);
  G4double amin = log(emin);
  G4double tgam, q;
  do {
    G4double x = amin + G4UniformRand()*(amax - amin);
    tgam = exp(x);
    q = amaj * G4UniformRand() / tgam;
  } while (q > a/tgam + b);

  // Sample gamma angle (Z - axis along the parent particle).
  // Universal distribution suggested by L. Urban (Geant3 manual (1993) 
  // Phys211) derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  G4double totalEnergy = kineticEnergy + electron_mass_c2;   

  const G4double a1 = 0.625, a2 = 3.*a1, d = 27.;
  G4double u = - log(G4UniformRand()*G4UniformRand());

  if (9./(9.+d) > G4UniformRand()) u /= a1;
  else                             u /= a2;
    
  G4double theta = u*electron_mass_c2/totalEnergy;
  G4double phi   = twopi * G4UniformRand();
  G4double dirz  = cos(theta);
  G4double sint  = sqrt(1. - dirz*dirz);
  G4double dirx  = sint*cos(phi);
  G4double diry  = sint*sin(phi); 
    
  G4ThreeVector gamDirection (dirx, diry, dirz);
  G4ThreeVector elecDirection = aParticle->GetMomentumDirection();
    
  gamDirection.rotateUz(elecDirection);   
  
  //
  // Update the incident particle 
  //
    
  G4double finalEnergy = kineticEnergy - tgam;  
    
  // Kinematic problem
  if (finalEnergy < 0.) {
    tgam += finalEnergy;
    finalEnergy = 0.0;
  }

  G4double mom = sqrt((totalEnergy + electron_mass_c2)*kineticEnergy);

  G4double finalX = mom*elecDirection.x() - tgam*gamDirection.x();
  G4double finalY = mom*elecDirection.y() - tgam*gamDirection.y();
  G4double finalZ = mom*elecDirection.z() - tgam*gamDirection.z();
      
  theChange->SetMomentumChange(finalX, finalY, finalZ);
  theChange->SetEnergyChange( finalEnergy );

  // create G4DynamicParticle object for the gamma 
  G4DynamicParticle* aGamma= new G4DynamicParticle (G4Gamma::Gamma(),
						    gamDirection, tgam);
	
  theChange->AddSecondary(aGamma); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LowEnergyBremsstrahlungGen::FindValueA(G4int Z, G4double e) const
{
  G4double value = 0.;
  
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> >::const_iterator pos;
  pos = paramA.find(Z);
  if (pos!= paramA.end())
    {
      G4VEMDataSet* dataSet = pos->second;
      value = dataSet->FindValue(e);
    }
  else
    {
      G4cout << "WARNING: G4LowEnergyBremsstrahlungGen::FindValue "
             << "did not find Z = "
	     << Z << G4endl;
    }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyBremsstrahlungGen::PrintData() const
{
  G4cout << generatorName 
         << ": e- Bremsstrahlung generation using EEDL database, "
         << "\n       gamma energy sampled from a parametrised formula."
         << "\n       Energy range = ["
         << lowEnergyLimit/eV << "eV, "
         << highEnergyLimit/GeV << "GeV]" 
         << G4endl;
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....












