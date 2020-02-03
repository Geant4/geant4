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
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created, with dummy implementation
// 12.09.01 V.Ivanchenko    Add param and interpolation of parameters
// 04.10.01 V.Ivanchenko    Add BindingEnergy method
// 25.10.01 MGP             Many bug fixes, mostly related to the
//                          management of pointers
// 29.11.01 V.Ivanchenko    New parametrisation + Excitation
// 30.05.02 V.Ivanchenko    Format and names of the data files were
//                          chenged to "ion-..."
// 17.02.04 V.Ivanchenko    Increase buffer size
// 23.03.09 L.Pandola       Updated warning message
// 03.12.10 V.Ivanchenko    Fixed memory leak in LoadData
//
// -------------------------------------------------------------------

#include "G4eIonisationParameters.hh"
#include "G4SystemOfUnits.hh"
#include "G4VEMDataSet.hh"
#include "G4ShellEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4LinInterpolation.hh"
#include "G4LogLogInterpolation.hh"
#include "G4LinLogLogInterpolation.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4Material.hh"
#include "G4DataVector.hh"
#include <fstream>
#include <sstream>


G4eIonisationParameters:: G4eIonisationParameters(G4int minZ, G4int maxZ)
  : zMin(minZ), zMax(maxZ),
  length(24)
{
  LoadData();
}


G4eIonisationParameters::~G4eIonisationParameters()
{ 
  // Reset the map of data sets: remove the data sets from the map 
  std::map<G4int,G4VEMDataSet*,std::less<G4int> >::iterator pos;

  for (pos = param.begin(); pos != param.end(); ++pos)
    {
      G4VEMDataSet* dataSet = (*pos).second;
      delete dataSet;
    }

  for (pos = excit.begin(); pos != excit.end(); ++pos)
    {
      G4VEMDataSet* dataSet = (*pos).second;
      delete dataSet;
    }

  activeZ.clear();
}


G4double G4eIonisationParameters::Parameter(G4int Z, G4int shellIndex, 
                                                     G4int parameterIndex, 
                                                     G4double e) const
{
  G4double value = 0.;
  G4int id = Z*100 + parameterIndex;
  std::map<G4int,G4VEMDataSet*,std::less<G4int> >::const_iterator pos;

  pos = param.find(id);
  if (pos!= param.end()) {
    G4VEMDataSet* dataSet = (*pos).second;
    G4int nShells = dataSet->NumberOfComponents();

    if(shellIndex < nShells) { 
      const G4VEMDataSet* component = dataSet->GetComponent(shellIndex);
      const G4DataVector ener = component->GetEnergies(0);
      G4double ee = std::max(ener.front(),std::min(ener.back(),e));
      value = component->FindValue(ee);
    } else {
      G4cout << "WARNING: G4IonisationParameters::FindParameter "
             << "has no parameters for shell= " << shellIndex 
             << "; Z= " << Z
             << G4endl;
    }
  } else {
    G4cout << "WARNING: G4IonisationParameters::Parameter "
           << "did not find ID = "
           << shellIndex << G4endl;
  }

  return value;
}

G4double G4eIonisationParameters::Excitation(G4int Z, G4double e) const
{
  G4double value = 0.;
  std::map<G4int,G4VEMDataSet*,std::less<G4int> >::const_iterator pos;

  pos = excit.find(Z);
  if (pos!= excit.end()) {
    G4VEMDataSet* dataSet = (*pos).second;

    const G4DataVector ener = dataSet->GetEnergies(0);
    G4double ee = std::max(ener.front(),std::min(ener.back(),e));
    value = dataSet->FindValue(ee);
  } else {
    G4cout << "WARNING: G4IonisationParameters::Excitation "
           << "did not find ID = "
           << Z << G4endl;
  }

  return value;
}


void G4eIonisationParameters::LoadData()
{
  // ---------------------------------------
  // Please document what are the parameters 
  // ---------------------------------------

  // define active elements

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  if (materialTable == 0)
     G4Exception("G4eIonisationParameters::LoadData",
		    "em1001",FatalException,"Unable to find MaterialTable");

  G4int nMaterials = G4Material::GetNumberOfMaterials();
  
  for (G4int mLocal=0; mLocal<nMaterials; mLocal++) {

    const G4Material* material= (*materialTable)[mLocal];        
    const G4ElementVector* elementVector = material->GetElementVector();
    const size_t nElements = material->GetNumberOfElements();
      
    for (size_t iEl=0; iEl<nElements; iEl++) {
      G4Element* element = (*elementVector)[iEl];
      G4double Z = element->GetZ();
      if (!(activeZ.contains(Z))) {
	activeZ.push_back(Z);
      }
    }
  }
  
  char* path = std::getenv("G4LEDATA");
  if (!path)
    { 
      G4Exception("G4BremsstrahlungParameters::LoadData",
		    "em0006",FatalException,"G4LEDATA environment variable not set");
      return;
    }
    
  G4String pathString(path);
  G4String path2("/ioni/ion-sp-");
  pathString += path2;  
  
  G4double energy, sum;
  
  size_t nZ = activeZ.size();
  
  for (size_t i=0; i<nZ; i++) {
    
    G4int Z = (G4int)activeZ[i];
    std::ostringstream ost;
    ost << pathString << Z << ".dat";
    G4String name(ost.str());

    std::ifstream file(name);
    std::filebuf* lsdp = file.rdbuf();

    if (! (lsdp->is_open()) ) {
      G4String excep = G4String("data file: ")
      + name + G4String(" not found");
      G4Exception("G4eIonisationParameters::LoadData",
		    "em0003",FatalException,excep);
    }

    // The file is organized into:
    // For each shell there are two lines:
    //    1st line:
    // 1st column is the energy of incident e-,
    // 2d  column is the parameter of screan term;
    //    2d line:
    // 3 energy (MeV) subdividing different approximation area of the spectrum
    // 20 point on the spectrum
    // The shell terminates with the pattern: -1   -1
    // The file terminates with the pattern: -2   -2

    std::vector<G4VEMDataSet*> p;
    for (size_t k=0; k<length; k++) 
      {
	G4VDataSetAlgorithm* inter = new G4LinLogLogInterpolation();
	G4VEMDataSet* composite = new G4CompositeEMDataSet(inter,1.,1.);
	p.push_back(composite);
      }

    G4int shell = 0;
    std::vector<G4DataVector*> a;
    a.resize(length);
    G4DataVector e;
    G4bool isReady = false;
    do {
      file >> energy >> sum;
      //Check if energy is valid
      if (energy < -2)
	{
	  G4String excep("invalid data file");
          G4Exception("G4eIonisationParameters::LoadData",
		    "em0005",FatalException,excep);
	  return;
	}

      if (energy == -2) break;

      if (energy >  -1) {
	// new energy
	if(!isReady) {
	  isReady = true;
	  e.clear();
	  for (size_t j=0; j<length; ++j) { 
	    a[j] = new G4DataVector(); 
	  }  
	}
        e.push_back(energy);
        a[0]->push_back(sum);
        for (size_t j=1; j<length; ++j) {
	  G4double qRead;
	  file >> qRead;
	  a[j]->push_back(qRead);
	}    

      } else {

        // End of set for a shell, fill the map
	for (size_t k=0; k<length; k++) {

          G4VDataSetAlgorithm* interp;
          if(0 == k) { interp  = new G4LinLogLogInterpolation(); }
          else       { interp  = new G4LogLogInterpolation(); }

	  G4DataVector* eVector = new G4DataVector;
	  size_t eSize = e.size();
	  for (size_t sLocal=0; sLocal<eSize; sLocal++) {
	    eVector->push_back(e[sLocal]);
	  }
	  G4VEMDataSet* set = new G4EMDataSet(shell,eVector,a[k],interp,1.,1.);

	  p[k]->AddComponent(set);
	} 
	  
	// next shell
        ++shell;
        isReady = false;
      }
    } while (energy > -2);
    
    file.close();

    for (size_t kk=0; kk<length; kk++) 
      {
	G4int id = Z*100 + kk;
	param[id] = p[kk];
      }
  }

  G4String pathString_a(path);
  G4String name_a = pathString_a + G4String("/ioni/ion-ex-av.dat");  
  std::ifstream file_a(name_a);
  std::filebuf* lsdp_a = file_a.rdbuf();
  G4String pathString_b(path);
  G4String name_b = pathString_b + G4String("/ioni/ion-ex-sig.dat");  
  std::ifstream file_b(name_b);
  std::filebuf* lsdp_b = file_b.rdbuf();
  
  if (! (lsdp_a->is_open()) ) {
     G4String excep = G4String("cannot open file ")
                    + name_a;
     G4Exception("G4eIonisationParameters::LoadData",
		    "em0003",FatalException,excep);
  }  
  if (! (lsdp_b->is_open()) ) {
     G4String excep = G4String("cannot open file ")
                    + name_b;
     G4Exception("G4eIonisationParameters::LoadData",
		    "em0003",FatalException,excep);
  }  

  // The file is organized into two columns:
  // 1st column is the energy
  // 2nd column is the corresponding value
  // The file terminates with the pattern: -1   -1
  //                                       -2   -2

  G4double ener, ener1, sig, sig1;
  G4int z    = 0;

  G4DataVector e;
  e.clear();
  G4DataVector d;
  d.clear();

  do {
    file_a >> ener >> sig;
    file_b >> ener1 >> sig1;
    
    if(ener != ener1) {
      G4cout << "G4eIonisationParameters: problem in excitation data "
             << "ener= " << ener
             << " ener1= " << ener1
             << G4endl;
    }

    // End of file
    if (ener == -2) {
      break;

      // End of next element
    } else if (ener == -1) {

      z++;
      G4double Z = (G4double)z;
    
	// fill map if Z is used
      if (activeZ.contains(Z)) {

	G4VDataSetAlgorithm* inter  = new G4LinInterpolation();
	G4DataVector* eVector = new G4DataVector;
	G4DataVector* dVector = new G4DataVector;
	size_t eSize = e.size();
	for (size_t sLocal2=0; sLocal2<eSize; sLocal2++) {
           eVector->push_back(e[sLocal2]);
           dVector->push_back(d[sLocal2]);
	}
	G4VEMDataSet* set = new G4EMDataSet(z,eVector,dVector,inter,1.,1.);
    	excit[z] = set;	
      }
      e.clear();
      d.clear();
  
    } else {

      e.push_back(ener*MeV);
      d.push_back(sig1*sig*barn*MeV);
    }
  } while (ener != -2);
  
  file_a.close();

}


void G4eIonisationParameters::PrintData() const
{
  G4cout << G4endl;
  G4cout << "===== G4eIonisationParameters =====" << G4endl;
  G4cout << G4endl;

  size_t nZ = activeZ.size();
  std::map<G4int,G4VEMDataSet*,std::less<G4int> >::const_iterator pos;

  for (size_t i=0; i<nZ; i++) {
    G4int Z = (G4int)activeZ[i];      

    for (size_t j=0; j<length; j++) {
    
      G4int index = Z*100 + j;

      pos = param.find(index);
      if (pos!= param.end()) {
        G4VEMDataSet* dataSet = (*pos).second;
        size_t nShells = dataSet->NumberOfComponents();

        for (size_t k=0; k<nShells; k++) {

          G4cout << "===== Z= " << Z << " shell= " << k 
                 << " parameter[" << j << "]  =====" 
                 << G4endl;
          const G4VEMDataSet* comp = dataSet->GetComponent(k);
          comp->PrintData();
	}
      }
    }
  }
  G4cout << "====================================" << G4endl;
}


