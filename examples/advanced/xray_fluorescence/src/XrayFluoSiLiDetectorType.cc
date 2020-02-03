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
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  19 Jun 2003  Alfonso Mantero Created
//
// -------------------------------------------------------------------

#include <fstream>
#include <sstream>

#include "XrayFluoSiLiDetectorType.hh"
#include "XrayFluoDataSet.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4DataVector.hh"
#include "G4LogLogInterpolation.hh"
#include "G4ios.hh"
#include "Randomize.hh"

XrayFluoSiLiDetectorType::XrayFluoSiLiDetectorType():
  detectorMaterial("SiLi"),efficiencySet(0)
{
  LoadResponseData("SILIresponse");

  LoadEfficiencyData("SILIefficiency");

 
}
XrayFluoSiLiDetectorType::~XrayFluoSiLiDetectorType()
{
  std::map<G4int,G4DataVector*,std::less<G4int> >::iterator pos;

  for (pos = energyMap.begin(); pos != energyMap.end(); pos++)
    {
      G4DataVector* dataSet = (*pos).second;
      delete dataSet;
      dataSet = 0;
    }
  for (pos = dataMap.begin(); pos != dataMap.end(); pos++)
    {
      G4DataVector* dataSet = (*pos).second;
      delete dataSet;
      dataSet = 0;
    }

  delete interpolation4; 

}
G4String XrayFluoSiLiDetectorType::GetDetectorMaterial()
{
  return detectorMaterial;
}

XrayFluoSiLiDetectorType* XrayFluoSiLiDetectorType::instance = 0;

XrayFluoSiLiDetectorType* XrayFluoSiLiDetectorType::GetInstance()

{
  if (instance == 0)
    {
      instance = new XrayFluoSiLiDetectorType;
     
    }
  return instance;
}


G4double XrayFluoSiLiDetectorType::ResponseFunction(G4double energy)
{
 
  G4double eMin = 1.500 *keV;
  G4double eMax = 6.403 *keV; 
  G4double value = 0.;
  G4double efficiency = 1.;
  
  const XrayFluoDataSet* dataSet = efficiencySet;
  G4int id = 0;
  G4DataVector energyVector;
  energyVector.push_back(1.486* keV);
  energyVector.push_back(1.740* keV);
  energyVector.push_back(3.688* keV);
  energyVector.push_back(4.510* keV);
  energyVector.push_back(5.414* keV);
  energyVector.push_back(6.404* keV);

  G4double infEnergy = 0 *keV;
  G4double supEnergy = 10* keV;
  G4int energyNumber = 0;
  G4double random = G4UniformRand();
 
  if (energy>=eMin && energy <=eMax)
    {

      for (G4int i=0; i<(G4int)energyVector.size(); i++){
	if (energyVector[i]/keV < energy/keV){

	  infEnergy = energyVector[i];
	  supEnergy = energyVector[i+1];

	  energyNumber = i+1;

	}
      }
      

      G4double infData = GetInfData(energy, random, energyNumber);

      G4double supData = GetSupData(energy,random, energyNumber);

      value = (std::log10(infData)*std::log10(supEnergy/energy) +
	       std::log10(supData)*std::log10(energy/infEnergy)) / 
	std::log10(supEnergy/infEnergy);
      value = std::pow(10,value);


    }
//    else if (energy<eMin || energy>eMax)
//      { 
//        G4double infEnergy = eMin;
//        G4double supEnergy = eMin/keV +1*keV;
 
//        G4double infData = GetInfData(eMin, random);
//        G4double supData = GetSupData(eMin,random);
//        value = (std::log10(infData)*std::log10(supEnergy/eMin) +
//  	       std::log10(supData)*std::log10(eMin/infEnergy)) / 
//  	std::log10(supEnergy/infEnergy);
//        value = std::pow(10,value);
//        value = value-eMin+ energy;


//      }


  else if (energy > eMax)
    {

      energyNumber = 5;

      value = (GetSupData(energy, random, energyNumber))+(energy - 6.404* keV);

    }



  else 
    {
 
      energyNumber = 1;

      value = (GetInfData(energy, random, energyNumber))+(energy - 1.486* keV);


      //      G4double mean = -14.03 * eV + 1.0047*energy/eV;
      //      G4double stdDev = 35.38 * eV + 0.004385*energy/eV;
      
      
      //      value = (G4RandGauss::shoot(mean,stdDev))*eV;

    }
  G4double  RandomNum = G4UniformRand(); 
  
  efficiency = dataSet->FindValue(value,id);
  if ( RandomNum>efficiency )
    {
      value = 0.;
    }

  //  G4cout << value << G4endl; 
  return value;

}
G4double XrayFluoSiLiDetectorType::GetInfData(G4double, G4double random, G4int posIndex)
{
  G4double value = 0.;
  G4int zMin = 1;
  G4int zMax = 6; 
  
  G4int Z = posIndex;
  
  if (Z<zMin) {Z=zMin;}
  if (Z>zMax) {Z=zMax;}
  
  if (Z >= zMin && Z <= zMax)
    {
      std::map<G4int,G4DataVector*,std::less<G4int> >::const_iterator pos;
      pos = energyMap.find(Z);
      std::map<G4int,G4DataVector*,std::less<G4int> >::const_iterator posData;
      posData = dataMap.find(Z);


      if (pos!= energyMap.end())
	{
	  G4DataVector energySet = *((*pos).second);
	  G4DataVector dataSet = *((*posData).second);
	  G4int nData = energySet.size();


	  G4double dataSum = 0;
	  G4double partSum = 0;
	  G4int index = 0;

	  // if data is not perfectly normalized (it may happen) 
	  // rnadom number is renormalized, in case it is higer 
	  //than the sum of all energies => segmentation fault.
	  
	  for  (G4int i = 0; i<nData; i++){
	    dataSum += dataSet[i];
	  }

	  G4double normRandom = random*dataSum;


	  while (normRandom> partSum)
	    {
	      
	      partSum += dataSet[index];
	      index++;
	    }
	 
	  
	  if (index >= 0 && index < nData)
	    {
	      value = energySet[index];
	      
	    }
	  
	}
    }
  return value;
}
G4double XrayFluoSiLiDetectorType::GetSupData(G4double, G4double random, G4int posIndex)
{
  G4double value = 0.;
  G4int zMin = 1;
  G4int zMax = 6;
  G4int Z = (posIndex+1);
  
  if (Z<zMin) {Z=zMin;}
  if (Z>zMax) {Z=zMax;}
  if (Z >= zMin && Z <= zMax)
    {
      std::map<G4int,G4DataVector*,std::less<G4int> >::const_iterator pos;
      pos = energyMap.find(Z);
      std::map<G4int,G4DataVector*,std::less<G4int> >::const_iterator posData;
      posData = dataMap.find(Z);
      if (pos!= energyMap.end())
	{
	  G4DataVector energySet = *((*pos).second);
	  G4DataVector dataSet = *((*posData).second);
	  G4int nData = energySet.size();
	  G4double dataSum = 0;
	  G4double partSum = 0;
	  G4int index = 0;

	  // if data is not perfectly normalized (it may happen) 
	  // rnadom number is renormalized, in case it is higer 
	  //than the sum of all energies => segmentation fault.

	  for  (G4int i = 0; i<nData; i++){
	    dataSum += dataSet[i];
	  }	  

	  G4double normRandom = random*dataSum;


	  while (normRandom> partSum)
	    {
	      partSum += dataSet[index];
	      index++;
	    }
	  
	  
	  if (index >= 0 && index < nData)
	    {
	      value = energySet[index];
	    }
	}
    }
  return value;
}
void XrayFluoSiLiDetectorType::LoadResponseData(G4String fileName)
{
  std::ostringstream ost;
  ost << fileName<<".dat";
  G4String name = ost.str();
  
  char* path = std::getenv("XRAYDATA");
  G4String dirFile;
  if (path) {
    G4String pathString(path);
    pathString += "\0";
    dirFile = pathString + "/" + name;
  }
  else{ 
    path = std::getenv("PWD");
    G4String pathString(path);
    pathString += "\0";
    dirFile = pathString + "/" + name;
  }
  
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
    {
      G4ExceptionDescription execp;
      execp <<  "XrayFluoSiLiDetectorType - data file: " + dirFile + " not found";
      G4Exception("XrayFluoSiLiDetectorType::LoadResponseData()","example-xray_fluorescence07",
	  FatalException, execp);
	}
  G4double a = 0;
  G4int k = 1;
  G4int q = 0;
  
  G4int Z = 1;
  G4DataVector* energies = new G4DataVector;
  G4DataVector* data = new G4DataVector;

  do
    {
      file >> a;
      G4int nColumns = 2;
      if (a == -1)
	{
	  if (q == 0)
	    {
	      // End of a  data set
	      energyMap[Z] = energies;
	      dataMap[Z] = data;
	      // Start of new shell data set
	      energies = new G4DataVector;
	      data = new G4DataVector;
	      Z++;	    
	    }      
	  q++;
	  if (q == nColumns)
	    {
	      q = 0;
	    }
	}
      else if (a == -2)
	{
	  // End of file; delete the empty vectors 
	  //created when encountering the last -1 -1 row
	  delete energies;
	  delete data;
	  
	}
      else
	{
	  // 1st column is energy
	  if(k%nColumns != 0)
	    {	
	      G4double e = a * keV;
	      energies->push_back(e);
	      k++;
	    }
	  else if (k%nColumns == 0)
	    {
	      // 2nd column is data
	      
	      data->push_back(a);
	      k = 1;
	    }
	}
    } while (a != -2); // end of file
  file.close();    
}

void XrayFluoSiLiDetectorType::LoadEfficiencyData(G4String fileName)
{
  char* path = std::getenv("XRAYDATA");
  G4String dirFile;
  if (path) {
    G4String pathString(path);
    dirFile = pathString + "/" + fileName;
  }
  else{ 
    path = std::getenv("PWD");
    G4String pathString(path);
    dirFile = pathString + "/" + fileName;
  }

  interpolation4 = new G4LogLogInterpolation();
  efficiencySet = new XrayFluoDataSet(1,fileName,interpolation4,keV,1);
}

