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
// History:
// -----------
// 31 Jul 2001   MGP        Created
//
// 15 Jul 2009   Nicolas A. Karakatsanis
//
//                           - New Constructor was added when logarithmic data are loaded as well
//                             to enhance CPU performance
//
//                           - Destructor was modified accordingly
//
//                           - LoadNonLogData method was created to load only the non-logarithmic data from G4EMLOW
//                             dataset. It is essentially performing the data loading operations as in the past.
//
//                           - LoadData method was revised in order to calculate the logarithmic values of the data
//                             It retrieves the data values from the G4EMLOW data files but, then, calculates the
//                             respective log values and loads them to seperate data structures. 
//
//                           - SetLogEnergiesData method was cretaed to set logarithmic values to G4 data vectors.
//                             The EM data sets, initialized this way, contain both non-log and log values.
//                             These initialized data sets can enhance the computing performance of data interpolation
//                             operations
//
// 26 Dec 2010 V.Ivanchenko Fixed Coverity warnings and cleanup logic
//
// -------------------------------------------------------------------

#include "G4EMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include <fstream>
#include <sstream>
#include "G4Integrator.hh"
#include "Randomize.hh"
#include "G4LinInterpolation.hh"

G4EMDataSet::G4EMDataSet(G4int Z, 
			 G4VDataSetAlgorithm* algo, 
			 G4double xUnit, 
			 G4double yUnit,
			 G4bool random):  
  energies(nullptr),
  data(nullptr),
  log_energies(nullptr),
  log_data(nullptr),
  algorithm(algo),
  pdf(nullptr),
  unitEnergies(xUnit),
  unitData(yUnit), 
  z(Z),
  randomSet(random)
{
  if (algorithm == nullptr) {
    G4Exception("G4EMDataSet::G4EMDataSet",
		    "em1012",FatalException,"interpolation == 0");
  } else if (randomSet) { BuildPdf(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EMDataSet::G4EMDataSet(G4int argZ, 
			 G4DataVector* dataX, 
			 G4DataVector* dataY, 
			 G4VDataSetAlgorithm* algo, 
			 G4double xUnit, 
			 G4double yUnit,
			 G4bool random):
  energies(dataX),
  data(dataY),
  log_energies(nullptr),
  log_data(nullptr),
  algorithm(algo),
  pdf(nullptr),
  unitEnergies(xUnit),
  unitData(yUnit), 
  z(argZ),
  randomSet(random)
{
  if (!algorithm || !energies || !data) {
    G4Exception("G4EMDataSet::G4EMDataSet",
		    "em1012",FatalException,"interpolation == 0");
  } else {
    if (energies->size() != data->size()) { 
      G4Exception("G4EMDataSet::G4EMDataSet",
		  "em1012",FatalException,"different size for energies and data");
    } else if (randomSet) {
      BuildPdf();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EMDataSet::G4EMDataSet(G4int argZ, 
			 G4DataVector* dataX, 
			 G4DataVector* dataY,
                         G4DataVector* dataLogX,
                         G4DataVector* dataLogY, 
			 G4VDataSetAlgorithm* algo, 
			 G4double xUnit, 
			 G4double yUnit,
			 G4bool random):
  energies(dataX),
  data(dataY),
  log_energies(dataLogX),
  log_data(dataLogY),
  algorithm(algo),
  pdf(nullptr),
  unitEnergies(xUnit),
  unitData(yUnit),
  z(argZ),
  randomSet(random)
{
  if (!algorithm || !energies || !data || !log_energies || !log_data) {
    G4Exception("G4EMDataSet::G4EMDataSet",
		    "em1012",FatalException,"interpolation == 0");
  } else {
    if ((energies->size() != data->size()) ||
	(energies->size() != log_energies->size()) ||
	(energies->size() != log_data->size())) { 
      G4Exception("G4EMDataSet::G4EMDataSet",
		  "em1012",FatalException,"different size for energies and data");
    } else if (randomSet) {
      BuildPdf();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EMDataSet::~G4EMDataSet()
{ 
  delete algorithm;
  delete energies; 
  delete data; 
  delete pdf; 
  delete log_energies; 
  delete log_data; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EMDataSet::FindValue(G4double energy, G4int /* componentId */) const
{
  if (energy <= (*energies)[0]) {
    return (*data)[0];
  }
  std::size_t i = energies->size()-1;
  if (energy >= (*energies)[i]) { return (*data)[i]; }

  //Nicolas A. Karakatsanis: Check if the logarithmic data have been loaded to decide
  //                         which Interpolation-Calculation method will be applied
  if (log_energies != nullptr) 
   {
     return algorithm->Calculate(energy, (G4int)FindLowerBound(energy), *energies, *data, *log_energies, *log_data);
   }

  return algorithm->Calculate(energy, (G4int)FindLowerBound(energy), *energies, *data);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EMDataSet::PrintData(void) const
{
  std::size_t size = energies->size();
  for (std::size_t i(0); i<size; ++i)
    {
      G4cout << "Point: " << ((*energies)[i]/unitEnergies)
	     << " - Data value: " << ((*data)[i]/unitData);
      if (pdf != 0) G4cout << " - PDF : " << (*pdf)[i];
      G4cout << G4endl; 
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EMDataSet::SetEnergiesData(G4DataVector* dataX, 
				  G4DataVector* dataY, 
				  G4int /* componentId */)
{
  if(!dataX || !dataY) {
    G4Exception("G4EMDataSet::SetEnergiesData",
		    "em1012",FatalException,"new interpolation == 0");
  } else {
    if (dataX->size() != dataY->size()) { 
      G4Exception("G4EMDataSet::SetEnergiesData",
		  "em1012",FatalException,"different size for energies and data");
    } else {

      delete energies; 
      energies = dataX;

      delete data; 
      data = dataY;
      //G4cout << "Size of energies: " << energies->size() << G4endl 
      //<< "Size of data: " << data->size() << G4endl;
    }
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EMDataSet::SetLogEnergiesData(G4DataVector* dataX,
                                     G4DataVector* dataY,
                                     G4DataVector* data_logX, 
				     G4DataVector* data_logY,
                                     G4int /* componentId */)
{
  if(!dataX || !dataY || !data_logX || !data_logY) {
    G4Exception("G4EMDataSet::SetEnergiesData",
		    "em1012",FatalException,"new interpolation == 0");
  } else {
    if (dataX->size() != dataY->size() ||
	dataX->size() != data_logX->size() ||
	dataX->size() != data_logY->size()) {
      G4Exception("G4EMDataSet::SetEnergiesData",
		  "em1012",FatalException,"different size for energies and data");
    } else {

      delete energies; 
      energies = dataX;

      delete data; 
      data = dataY;

      delete log_energies; 
      log_energies = data_logX;

      delete log_data; 
      log_data = data_logY;
      //G4cout << "Size of energies: " << energies->size() << G4endl 
      //<< "Size of data: " << data->size() << G4endl;
    }
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EMDataSet::LoadData(const G4String& fileName)
{
 // The file is organized into four columns:
 // 1st column contains the values of energy
 // 2nd column contains the corresponding data value
 // The file terminates with the pattern: -1   -1
 //                                       -2   -2

  G4String fullFileName(FullFileName(fileName));
  std::ifstream in(fullFileName);

  if (!in.is_open())
    {
      G4String message("data file \"");
      message += fullFileName;
      message += "\" not found";
      G4Exception("G4EMDataSet::LoadData",
		    "em1012",FatalException,message);
      return false;
    }

  delete energies;
  delete data;   
  delete log_energies;
  delete log_data;   
  energies = new G4DataVector;
  data = new G4DataVector;
  log_energies = new G4DataVector;
  log_data = new G4DataVector;

  G4double a, b;
  do
    {
      in >> a >> b;
  
      if (a != -1 && a != -2)
	{
	  if (a==0.) { a=1e-300; }
	  if (b==0.) { b=1e-300; }
          a *= unitEnergies;
          b *= unitData;
	  energies->push_back(a);
	  log_energies->push_back(std::log10(a));
	  data->push_back(b);
	  log_data->push_back(std::log10(b));
	}
    }
  while (a != -2);

  if (randomSet) { BuildPdf(); }
 
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EMDataSet::LoadNonLogData(const G4String& fileName)
{
 // The file is organized into four columns:
 // 1st column contains the values of energy
 // 2nd column contains the corresponding data value
 // The file terminates with the pattern: -1   -1
 //                                       -2   -2

  G4String fullFileName(FullFileName(fileName));
  std::ifstream in(fullFileName);
  if (!in.is_open())
    {
      G4String message("data file \"");
      message += fullFileName;
      message += "\" not found";
      G4Exception("G4EMDataSet::LoadNonLogData",
		    "em1012",FatalException,message);
    }

  G4DataVector* argEnergies=new G4DataVector;
  G4DataVector* argData=new G4DataVector;

  G4double a;
  G4int k = 0;
  G4int nColumns = 2;

  do
    {
      in >> a;
  
      if (a != -1 && a != -2)
	{
	  if (k%nColumns == 0)
            {
	     argEnergies->push_back(a*unitEnergies);
            }
	  else if (k%nColumns == 1)
            {
	     argData->push_back(a*unitData);
            }
          k++;
	}
    }
  while (a != -2);
 
  SetEnergiesData(argEnergies, argData, 0);

  if (randomSet) BuildPdf();
 
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EMDataSet::SaveData(const G4String& name) const
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
      G4String message("cannot open \"");
      message+=fullFileName;
      message+="\"";
      G4Exception("G4EMDataSet::SaveData",
		    "em1012",FatalException,message);
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::size_t G4EMDataSet::FindLowerBound(G4double x) const
{
  std::size_t lowerBound = 0;
  std::size_t upperBound(energies->size() - 1);
  
  while (lowerBound <= upperBound) 
    {
      std::size_t midBin((lowerBound + upperBound) / 2);

      if (x < (*energies)[midBin]) upperBound = midBin - 1;
      else lowerBound = midBin + 1;
    }
  
  return upperBound;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::size_t G4EMDataSet::FindLowerBound(G4double x, G4DataVector* values) const
{
  std::size_t lowerBound = 0;;
  std::size_t upperBound(values->size() - 1);
  
  while (lowerBound <= upperBound) 
    {
      std::size_t midBin((lowerBound + upperBound) / 2);

      if (x < (*values)[midBin]) upperBound = midBin - 1;
      else lowerBound = midBin + 1;
    }
  
  return upperBound;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String G4EMDataSet::FullFileName(const G4String& name) const
{
  const char* path = G4FindDataDir("G4LEDATA");
  if (!path) {
     G4Exception("G4EMDataSet::FullFileName",
		    "em0006",FatalException,"G4LEDATA environment variable not set");
    return "";
  }
  std::ostringstream fullFileName;
  fullFileName << path << '/' << name << z << ".dat";
                      
  return G4String(fullFileName.str().c_str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EMDataSet::BuildPdf() 
{
  pdf = new G4DataVector;
  G4Integrator <G4EMDataSet, G4double(G4EMDataSet::*)(G4double)> integrator;

  std::size_t nData = data->size();
  pdf->push_back(0.);

  // Integrate the data distribution 
  std::size_t i;
  G4double totalSum = 0.;
  for (i=1; i<nData; ++i)
    {
      G4double xLow = (*energies)[i-1];
      G4double xHigh = (*energies)[i];
      G4double sum = integrator.Legendre96(this, &G4EMDataSet::IntegrationFunction, xLow, xHigh);
      totalSum = totalSum + sum;
      pdf->push_back(totalSum);
    }

  // Normalize to the last bin
  G4double tot = 0.;
  if (totalSum > 0.) tot = 1. / totalSum;
  for (i=1;  i<nData; ++i)
    {
      (*pdf)[i] = (*pdf)[i] * tot;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EMDataSet::RandomSelect(G4int /* componentId */) const
{
  G4double value = 0.;
  // Random select a X value according to the cumulative probability distribution
  // derived from the data

  if (!pdf) {
    G4Exception("G4EMDataSet::RandomSelect",
		    "em1012",FatalException,"PDF has not been created for this data set");
    return value;
  }

  G4double x = G4UniformRand();

  // Locate the random value in the X vector based on the PDF
  G4int bin = (G4int)FindLowerBound(x,pdf);

  // Interpolate the PDF to calculate the X value: 
  // linear interpolation in the first bin (to avoid problem with 0),
  // interpolation with associated data set algorithm in other bins
  G4LinInterpolation linearAlgo;
  if (bin == 0) value = linearAlgo.Calculate(x, bin, *pdf, *energies);
  else value = algorithm->Calculate(x, bin, *pdf, *energies);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EMDataSet::IntegrationFunction(G4double x)
{
  // This function is needed by G4Integrator to calculate the integral of the data distribution
  G4double y = 0;

  // Locate the random value in the X vector based on the PDF
  G4int bin = (G4int)FindLowerBound(x);

  // Interpolate to calculate the X value: 
  // linear interpolation in the first bin (to avoid problem with 0),
  // interpolation with associated algorithm in other bins

  G4LinInterpolation linearAlgo;
  
  if (bin == 0) y = linearAlgo.Calculate(x, bin, *energies, *data);
  else y = algorithm->Calculate(x, bin, *energies, *data);
 
  return y;
}

