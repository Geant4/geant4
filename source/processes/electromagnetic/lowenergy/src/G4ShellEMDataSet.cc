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
//  1 Aug 2001   MGP         Created
//
//  09.10.01   V.Ivanchenko  Add case z=0
//
//  9 Mar 2008   MGP         Cleaned up unreadable code modified by former developer
//                           (Further clean-up needed)
//
//  15 Jul 2009   Nicolas A. Karakatsanis
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
//
// -------------------------------------------------------------------

#include "G4ShellEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include <fstream>
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4ShellEMDataSet::G4ShellEMDataSet(G4int zeta, G4VDataSetAlgorithm* algo, 
				   G4double eUnit, 
				   G4double dataUnit)
  : 
  algorithm(algo),
  unitEnergies(eUnit),
  unitData(dataUnit),
  z(zeta)
{
  if (algorithm == nullptr) 
    G4Exception("G4ShellEMDataSet::G4ShellEMDataSet()","em0007",
		FatalErrorInArgument, "Interpolation == 0");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4ShellEMDataSet::~G4ShellEMDataSet()
{
  CleanUpComponents();
  if (algorithm) delete algorithm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4ShellEMDataSet::FindValue(G4double energy, G4int /* componentId */) const
{
  // Returns the sum over the shells corresponding to e
  G4double value = 0.;

  std::vector<G4VEMDataSet *>::const_iterator i(components.begin());
  std::vector<G4VEMDataSet *>::const_iterator end(components.end());

  while (i != end)
    {
      value += (*i)->FindValue(energy);
      i++;
    }

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ShellEMDataSet::PrintData() const
{
  const G4int n = (G4int)NumberOfComponents();

  G4cout << "The data set has " << n << " components" << G4endl;
  G4cout << G4endl;
 
  G4int i = 0;
 
  while (i < n)
    {
      G4cout << "--- Component " << i << " ---" << G4endl;
      GetComponent(i)->PrintData();
      ++i;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ShellEMDataSet::SetEnergiesData(G4DataVector* energies, 
				       G4DataVector* data, 
				       G4int componentId)
{
  G4VEMDataSet* component = components[componentId]; 
  if (component)
    {
      component->SetEnergiesData(energies, data, 0);
      return;
    }

  G4String msg = "component ";
  msg += static_cast<char>(componentId);
  msg += " not found";
 
  G4Exception("G4ShellEMDataSet::SetEnergiesData()","em0008", FatalErrorInArgument ,msg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ShellEMDataSet::SetLogEnergiesData(G4DataVector* energies,
                                          G4DataVector* data,
                                          G4DataVector* log_energies,
                                          G4DataVector* log_data,
                                          G4int componentId)
{
  G4VEMDataSet* component = components[componentId]; 
  if (component)
    {
      component->SetLogEnergiesData(energies, data, log_energies, log_data, 0);
      return;
    }

  G4String msg = "component ";
  msg += static_cast<char>(componentId);
  msg += " not found";
  G4Exception("G4ShellEMDataSet::SetLogEnergiesData()","em0008", 
	      FatalErrorInArgument ,msg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4ShellEMDataSet::LoadData(const G4String& file)
{
  CleanUpComponents();

  G4String fullFileName = FullFileName(file);
  std::ifstream in(fullFileName);

  if (!in.is_open())
    {
      G4String message("Data file \"");
      message += fullFileName;
      message += "\" not found";
      G4Exception("G4ShellEMDataSet::LoadData()", "em0003",FatalException, message);
      return 0;
    }
  
  G4DataVector* orig_shell_energies = 0;
  G4DataVector* orig_shell_data = 0;
  G4DataVector* log_shell_energies = 0;
  G4DataVector* log_shell_data = 0;

  G4double a = 0.;
  G4int shellIndex = 0;
  G4int k = 0;
  G4int nColumns = 2;

  do
    {
      in >> a;

      if (a==0.) a=1e-300;

      // The file is organized into four columns:
      // 1st column contains the values of energy
      // 2nd column contains the corresponding data value
      // The file terminates with the pattern: -1   -1
      //                                       -2   -2
      //
      if (a == -1)
	{
	  if ((k%nColumns == 0) && (orig_shell_energies != 0) )
	    {
	     AddComponent(new G4EMDataSet(shellIndex, orig_shell_energies, 
					  orig_shell_data, log_shell_energies, 
					  log_shell_data, algorithm->Clone(), 
					  unitEnergies, unitData));
	      orig_shell_energies = 0;
	      orig_shell_data = 0;
              log_shell_energies = 0;
              log_shell_data = 0;
	    }
	}
      else if (a != -2)
	{
	  if (orig_shell_energies == 0)
	    {
	     orig_shell_energies = new G4DataVector;
	     orig_shell_data = new G4DataVector;
             log_shell_energies = new G4DataVector;
             log_shell_data = new G4DataVector;
	    }
	  if (k%nColumns == 0)
            {
	     orig_shell_energies->push_back(a*unitEnergies);
	     log_shell_energies->push_back(std::log10(a) + std::log10(unitEnergies));
            }
	  else if (k%nColumns == 1)
            {
	     orig_shell_data->push_back(a*unitData);
             log_shell_data->push_back(std::log10(a) + std::log10(unitData));
            }
          k++;
	}
      else k = 1;
    }
  while (a != -2);  // End of file

  delete orig_shell_energies;
  delete orig_shell_data;
  delete log_shell_energies;
  delete log_shell_data;

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4ShellEMDataSet::LoadNonLogData(const G4String& file)
{
  CleanUpComponents();

  G4String fullFileName = FullFileName(file);
  std::ifstream in(fullFileName);

  if (!in.is_open())
    {
      G4String message("G4ShellEMDataSet::LoadData - data file \"");
      message += fullFileName;
      message += "\" not found";
      G4Exception("G4ShellEMDataSet::LoadNonLogData()", "em0003",FatalException, message);
      return 0;
    }

  G4DataVector* orig_shell_energies = 0;
  G4DataVector* orig_shell_data = 0;

  G4double a = 0.;
  G4int shellIndex = 0;
  G4int k = 0;
  G4int nColumns = 2;

  do
    {
      in >> a;
      // The file is organized into four columns:
      // 1st column contains the values of energy
      // 2nd column contains the corresponding data value
      // The file terminates with the pattern: -1   -1
      //                                       -2   -2
      //
      if (a == -1)
	{
	  if ((k%nColumns == 0) && (orig_shell_energies != 0) )
	    {
	     AddComponent(new G4EMDataSet(shellIndex, orig_shell_energies, 
					  orig_shell_data, algorithm->Clone(), 
					  unitEnergies, unitData));
	      orig_shell_energies = 0;
	      orig_shell_data = 0;
	    }
	}
      else if (a != -2)
	{
	  if (orig_shell_energies == 0)
	    {
	     orig_shell_energies = new G4DataVector;
	     orig_shell_data = new G4DataVector;
	    }
	  if (k%nColumns == 0)
            {
	     orig_shell_energies->push_back(a*unitEnergies);
            }
	  else if (k%nColumns == 1)
            {
	     orig_shell_data->push_back(a*unitData);
            }
          k++;
	}
      else k = 1;
    }
  while (a != -2);  // End of file


  delete orig_shell_energies;
  delete orig_shell_data;

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4ShellEMDataSet::SaveData(const G4String& file) const
{
  G4String fullFileName = FullFileName(file);
  std::ofstream out(fullFileName);

  if (!out.is_open())
    {
      G4String message("Cannot open \"");
      message += fullFileName;
      message += "\"";
      G4Exception("G4EMDataSet::SaveData()","em0005",FatalException,message);
    }
 
  const G4int n = (G4int)NumberOfComponents();
  G4int k = 0;
 
  while (k < n)
    {
      const G4VEMDataSet* component = GetComponent(k);
  
      if (component)
	{
	  const G4DataVector& energies = component->GetEnergies(0);
	  const G4DataVector& data = component->GetData(0);
 
	  G4DataVector::const_iterator i = energies.cbegin();
	  G4DataVector::const_iterator endI = energies.cend();
	  G4DataVector::const_iterator j = data.cbegin();
  
	  while (i != endI)
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
  
      k++;
    }
 
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

void G4ShellEMDataSet::CleanUpComponents()
{
  while (!components.empty())
    {
      if (components.back()) delete components.back();
      components.pop_back();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String G4ShellEMDataSet::FullFileName(const G4String& fileName) const
{
  const char* path = G4FindDataDir("G4LEDATA");
  
  if (!path)
  {
    G4Exception("G4ShellEMDataSet::FullFileName()","em0006",JustWarning,"Please set G4LEDATA");
    return "";
  }
  
  std::ostringstream fullFileName;
 
  fullFileName << path << '/' << fileName << z << ".dat";
                      
  return G4String(fullFileName.str().c_str());
}
