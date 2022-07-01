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
// 31 Jul 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "G4RDEMDataSet.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include <fstream>
#include <sstream>
#include "G4Integrator.hh"
#include "Randomize.hh"
#include "G4RDLinInterpolation.hh"


G4RDEMDataSet::G4RDEMDataSet(G4int Z, 
			 G4RDVDataSetAlgorithm* algo, 
			 G4double xUnit, 
			 G4double yUnit,
			 G4bool random): 
  z(Z),
  energies(0),
  data(0),
  algorithm(algo),
  unitEnergies(xUnit),
  unitData(yUnit),
  pdf(0),
  randomSet(random)
{
  if (algorithm == 0)
    G4Exception("G4RDEMDataSet::G4RDEMDataSet()",
                "InvalidSetup", FatalException, "Interpolation == 0!");
  if (randomSet) BuildPdf();
}

G4RDEMDataSet::G4RDEMDataSet(G4int argZ, 
			 G4DataVector* dataX, 
			 G4DataVector* dataY, 
			 G4RDVDataSetAlgorithm* algo, 
			 G4double xUnit, 
			 G4double yUnit,
			 G4bool random):
  z(argZ),
  energies(dataX),
  data(dataY),
  algorithm(algo),
  unitEnergies(xUnit),
  unitData(yUnit),
  pdf(0),
  randomSet(random)
{
  if (algorithm == 0)
    G4Exception("G4RDEMDataSet::G4RDEMDataSet()",
                "InvalidSetup", FatalException, "Interpolation == 0!");

  if ((energies == 0) ^ (data == 0))
    G4Exception("G4RDEMDataSet::G4RDEMDataSet()", "InvalidSetup",
           FatalException, "Different size for energies and data (zero case)!");

  if (energies == 0) return;
  
  if (energies->size() != data->size()) 
    G4Exception("G4RDEMDataSet::G4RDEMDataSet()", "InvalidSetup",
           FatalException, "Different size for energies and data!");

  if (randomSet) BuildPdf();
}

G4RDEMDataSet::~G4RDEMDataSet()
{ 
  delete algorithm;
  if (energies) delete energies;
  if (data) delete data;
  if (pdf) delete pdf;
}

G4double G4RDEMDataSet::FindValue(G4double energy, G4int /* componentId */) const
{
  if (!energies)
    G4Exception("G4RDEMDataSet::FindValue()", "InvalidSetup",
                FatalException, "Energies == 0!");
  if (energies->empty()) return 0;
  if (energy <= (*energies)[0]) return (*data)[0];

  size_t i = energies->size()-1;
  if (energy >= (*energies)[i]) return (*data)[i];

  return algorithm->Calculate(energy, FindLowerBound(energy), *energies, *data);
}


void G4RDEMDataSet::PrintData(void) const
{
  if (!energies)
    {
      G4cout << "Data not available." << G4endl;
    }
  else
    {
      size_t size = energies->size();
      for (size_t i(0); i<size; i++)
	{
	  G4cout << "Point: " << ((*energies)[i]/unitEnergies)
		 << " - Data value: " << ((*data)[i]/unitData);
	  if (pdf != 0) G4cout << " - PDF : " << (*pdf)[i];
	  G4cout << G4endl; 
	}
    }
}


void G4RDEMDataSet::SetEnergiesData(G4DataVector* dataX, 
				  G4DataVector* dataY, 
				  G4int /* componentId */)
{
  if (energies) delete energies;
  energies = dataX;

  if (data) delete data;
  data = dataY;
 
  if ((energies == 0) ^ (data==0)) 
    G4Exception("G4RDEMDataSet::SetEnergiesData()", "InvalidSetup",
       FatalException, "Different size for energies and data (zero case)!");

  if (energies == 0) return;
  
  if (energies->size() != data->size()) 
    G4Exception("G4RDEMDataSet::SetEnergiesData()", "InvalidSetup",
       FatalException, "Different size for energies and data!");
}

G4bool G4RDEMDataSet::LoadData(const G4String& fileName)
{
  // The file is organized into two columns:
  // 1st column is the energy
  // 2nd column is the corresponding value
  // The file terminates with the pattern: -1   -1
  //                                       -2   -2
 
  G4String fullFileName(FullFileName(fileName));
  std::ifstream in(fullFileName);

  if (!in.is_open())
    {
      G4String message("Data file \"");
      message += fullFileName;
      message += "\" not found";
      G4Exception("G4RDEMDataSet::LoadData()", "DataNotFound",
                  FatalException, message);
    }

  G4DataVector* argEnergies=new G4DataVector;
  G4DataVector* argData=new G4DataVector;

  G4double a;
  bool energyColumn(true);

  do
    {
      in >> a;
  
      if (a!=-1 && a!=-2)
	{
	  if (energyColumn)
	    argEnergies->push_back(a*unitEnergies);
	  else
	    argData->push_back(a*unitData);
	  energyColumn=(!energyColumn);
	}
    }
  while (a != -2);
 
  SetEnergiesData(argEnergies, argData, 0);
  if (randomSet) BuildPdf();
 
  return true;
}

G4bool G4RDEMDataSet::SaveData(const G4String& name) const
{
  // The file is organized into two columns:
  // 1st column is the energy
  // 2nd column is the corresponding value
  // The file terminates with the pattern: -1   -1
  //                                       -2   -2
 
  G4String fullFileName(FullFileName(name));
  std::ofstream out(fullFileName);

  if (!out.is_open())
    {
      G4String message("Cannot open \"");
      message+=fullFileName;
      message+="\"";
      G4Exception("G4RDEMDataSet::SaveData()", "CannotOpenFile",
                  FatalException, message);
    }
 
  out.precision(10);
  out.width(15);
  out.setf(std::ofstream::left);
  
  if (energies!=0 && data!=0)
    {
      G4DataVector::const_iterator i(energies->begin());
      G4DataVector::const_iterator endI(energies->end());
      G4DataVector::const_iterator j(data->begin());
  
      while (i!=endI)
	{
	  out.precision(10);
	  out.width(15);
	  out.setf(std::ofstream::left);
	  out << ((*i)/unitEnergies) << ' ';

	  out.precision(10);
	  out.width(15);
	  out.setf(std::ofstream::left);
	  out << ((*j)/unitData) << std::endl;

	  i++;
	  j++;
	}
    }
 
  out.precision(10);
  out.width(15);
  out.setf(std::ofstream::left);
  out << -1.f << ' ';

  out.precision(10);
  out.width(15);
  out.setf(std::ofstream::left);
  out << -1.f << std::endl;

  out.precision(10);
  out.width(15);
  out.setf(std::ofstream::left);
  out << -2.f << ' ';

  out.precision(10);
  out.width(15);
  out.setf(std::ofstream::left);
  out << -2.f << std::endl;

  return true;
}

size_t G4RDEMDataSet::FindLowerBound(G4double x) const
{
  size_t lowerBound = 0;
  size_t upperBound(energies->size() - 1);
  
  while (lowerBound <= upperBound) 
    {
      size_t midBin((lowerBound + upperBound) / 2);

      if (x < (*energies)[midBin]) upperBound = midBin - 1;
      else lowerBound = midBin + 1;
    }
  
  return upperBound;
}


size_t G4RDEMDataSet::FindLowerBound(G4double x, G4DataVector* values) const
{
  size_t lowerBound = 0;;
  size_t upperBound(values->size() - 1);
  
  while (lowerBound <= upperBound) 
    {
      size_t midBin((lowerBound + upperBound) / 2);

      if (x < (*values)[midBin]) upperBound = midBin - 1;
      else lowerBound = midBin + 1;
    }
  
  return upperBound;
}


G4String G4RDEMDataSet::FullFileName(const G4String& name) const
{
  const char* path = G4FindDataDir("G4LEDATA");
  if (!path)
    G4Exception("G4RDEMDataSet::FullFileName()", "InvalidSetup",
                FatalException, "G4LEDATA environment variable not set!");
  
  std::ostringstream fullFileName;
  fullFileName << path << '/' << name << z << ".dat";
                      
  return G4String(fullFileName.str().c_str());
}


void G4RDEMDataSet::BuildPdf() 
{
  pdf = new G4DataVector;
  G4Integrator <G4RDEMDataSet, G4double(G4RDEMDataSet::*)(G4double)> integrator;

  G4int nData = data->size();
  pdf->push_back(0.);

  // Integrate the data distribution 
  G4int i;
  G4double totalSum = 0.;
  for (i=1; i<nData; i++)
    {
      G4double xLow = (*energies)[i-1];
      G4double xHigh = (*energies)[i];
      G4double sum = integrator.Legendre96(this, &G4RDEMDataSet::IntegrationFunction, xLow, xHigh);
      totalSum = totalSum + sum;
      pdf->push_back(totalSum);
    }

  // Normalize to the last bin
  G4double tot = 0.;
  if (totalSum > 0.) tot = 1. / totalSum;
  for (i=1;  i<nData; i++)
    {
      (*pdf)[i] = (*pdf)[i] * tot;
    }
}


G4double G4RDEMDataSet::RandomSelect(G4int /* componentId */) const
{
  // Random select a X value according to the cumulative probability distribution
  // derived from the data

  if (!pdf) G4Exception("G4RDEMDataSet::RandomSelect()", "InvalidSetup",
           FatalException, "PDF has not been created for this data set");

  G4double value = 0.;
  G4double x = G4UniformRand();

  // Locate the random value in the X vector based on the PDF
  size_t bin = FindLowerBound(x,pdf);

  // Interpolate the PDF to calculate the X value: 
  // linear interpolation in the first bin (to avoid problem with 0),
  // interpolation with associated data set algorithm in other bins

  G4RDLinInterpolation linearAlgo;
  if (bin == 0) value = linearAlgo.Calculate(x, bin, *pdf, *energies);
  else value = algorithm->Calculate(x, bin, *pdf, *energies);

  //  G4cout << x << " random bin "<< bin << " - " << value << G4endl;
  return value;
}

G4double G4RDEMDataSet::IntegrationFunction(G4double x)
{
  // This function is needed by G4Integrator to calculate the integral of the data distribution

  G4double y = 0;

 // Locate the random value in the X vector based on the PDF
  size_t bin = FindLowerBound(x);

  // Interpolate to calculate the X value: 
  // linear interpolation in the first bin (to avoid problem with 0),
  // interpolation with associated algorithm in other bins

  G4RDLinInterpolation linearAlgo;
  
  if (bin == 0) y = linearAlgo.Calculate(x, bin, *energies, *data);
  else y = algorithm->Calculate(x, bin, *energies, *data);
 
  return y;
}

