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
// $Id: G4CrossSectionChargeTransferExp.cc,v 1.6 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// Reference: TNS Geant4-DNA paper
// Reference for implementation model: NIM. 155, pp. 145-156, 1978

// History:
// -----------
// Date         Name              Modification
// 28 Dec 2007  M.G. Pia          Created
//
// -------------------------------------------------------------------

// Class description:
// Total cross section for incident p charge transfer from experimental data 
// interpolation
// Reference: http://www-cfadc.phy.ornl.gov/astro/ps/data/
//
// The fit formula is:
// ln(cross section)=\sum C_i T_i(x)
//            x=[2ln(E)-ln(E_min)-ln(E_max)]/[ln(E_max)-ln(E_min)]
//            T_1(x)=1
//            T_2(x)=x
//          T_n+2(x)=2T_n+1(x)-T_n(x)
//        
// Where cross section is given in 10^-16 cm2 and E is in keV,
// E_min=1e-4 keV and E_max=1.e3 keV.

// Further documentation available from http://www.ge.infn.it/geant4/

// -------------------------------------------------------------------

#include <sstream>

#include "G4CrossSectionChargeTransferExp.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4Material.hh"
#include "G4DataVector.hh"

G4CrossSectionChargeTransferExp::G4CrossSectionChargeTransferExp()
{
  // Default energy limits (defined for protection against anomalous behaviour only)
  name = "ChargeTransferExp";
  lowEnergyLimit = 0.1 * eV;
  highEnergyLimit = 1.0 * MeV;

  // Units: energies in eV, cross sections in cm2
  unit1 = eV;
  unit2 = cm2;  

  // Create DataSets from experimental data to be interpolated
  G4VEMDataSet* dataSet = 0;
  G4String materialName;

  // ---- MGP ---- Primitive implementation in 1st iteration
  // Load all experimental data sets available
  // Improvement in later iterations: load only the data sets of materials present in the
  // MaterialTable (i.e. actually used in the simulation), or load on demand

  // H2 cross section data
  dataSet = LoadData("/charge_transf/cs-p-H2");
  materialName = "H2";
  crossMap[materialName] = dataSet;

  // He cross section data
  dataSet = LoadData("/charge_transf/cs-p-He");
  materialName = "He";
  crossMap[materialName] = dataSet;

  // CO cross section data
  dataSet = LoadData("/charge_transf/cs-p-CO");
  materialName = "CO";
  crossMap[materialName] = dataSet;

  // CO2 cross section data
  dataSet = LoadData("/charge_transf/cs-p-CO2");
  materialName = "CO2";
  crossMap[materialName] = dataSet;

  // N2 cross section data
  dataSet = LoadData("/charge_transf/cs-p-N2");
  materialName = "N2";
  crossMap[materialName] = dataSet;
}


G4CrossSectionChargeTransferExp::~G4CrossSectionChargeTransferExp()
{
  std::map<G4String,G4VEMDataSet*,std::less<G4String> >::iterator pos;
  for (pos = crossMap.begin(); pos != crossMap.end(); ++pos)
    {
      G4VEMDataSet* dataSet = (*pos).second;
      delete dataSet;
      dataSet = 0;
    }
}


G4double G4CrossSectionChargeTransferExp::CrossSection(const G4Track& track)
{
  G4double sigma = DBL_MIN;
 
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double e = particle->GetKineticEnergy();

  // Get the material the track is in
  G4Material* material = track.GetMaterial();
  G4String materialName;
  materialName = material->GetName();

  // Check whether cross section data are available for the current material

  std::map<G4String,G4VEMDataSet*,std::less<G4String> >::const_iterator pos;
  pos = crossMap.find(materialName);
  if (pos!= crossMap.end())
    {
      G4VEMDataSet* dataSet = (*pos).second;
      // G4VEMDataSet* dataSet = crossMap["He"];
      const G4DataVector energies = dataSet->GetEnergies(0);
      G4int nPoints = energies.size(); 

      G4double eMin = energies[0];
      G4double eMax = energies[nPoints-1];

      if (e >=  eMin && e <= eMax)
	{
	  sigma = dataSet->FindValue(e);
	}      
    }
  return sigma;
}


G4VEMDataSet* G4CrossSectionChargeTransferExp::LoadData(const G4String& fileName)
{
  std::ostringstream ost;
  ost << fileName << ".dat";
  G4String nameF(ost.str());

  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep("G4EMDataSet - G4LEDATA environment variable not set");
      G4Exception(excep);
    }
  
  G4String pathString(path);
  G4String dirFile = pathString + nameF;
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();

  if (! (lsdp->is_open()) )
    {
      G4String s1("G4CrossSectionChargeTransferExp::LoadData data file: ");
      G4String s2(" not found");
      G4String excep = s1 + dirFile + s2;
      G4Exception(excep);
    }

  G4double a = 0;
  G4int k = 1;
  G4DataVector* energies = new G4DataVector;
  G4DataVector* data = new G4DataVector;
  do
    {
      file >> a;
      G4int nColumns = 2;
      // The file is organized into two columns:
      // 1st column is the energy
      // 2nd column is the corresponding value
      // The file terminates with the pattern: -1   -1
      //                                       -2   -2
      if (a == -1 || a == -2)
	{
	}
      else
	{
	  if (k%nColumns != 0)
	    {	
	      G4double e = a * unit1;
	      energies->push_back(e);
	      k++;
	    }
	  else if (k%nColumns == 0)
	    {
	      G4double value = a * unit2;
	      data->push_back(value);
	      k = 1;
	    }
	}
    } while (a != -2); // end of file
  
  file.close();
  G4VDataSetAlgorithm* algo = new G4LogLogInterpolation;
  G4VEMDataSet* dataSet = new G4EMDataSet(0,energies,data,algo);

  //dataSet->PrintData();

  return dataSet;

}
