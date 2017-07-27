
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
//  MR - 04/04/2012
//  Based on G4DNACrossSectionDataSet
// 


#include "G4MicroElecCrossSectionDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4EMDataSet.hh"
#include <vector>
#include <fstream>
#include <sstream>


G4MicroElecCrossSectionDataSet::G4MicroElecCrossSectionDataSet(G4VDataSetAlgorithm* argAlgorithm, 
						   G4double argUnitEnergies, 
						   G4double argUnitData)
  :
   algorithm(argAlgorithm), unitEnergies(argUnitEnergies), unitData(argUnitData)
{
  z = 0;

}


G4MicroElecCrossSectionDataSet::~G4MicroElecCrossSectionDataSet()
{
  CleanUpComponents();
 
  if (algorithm)
    delete algorithm;
}

G4bool G4MicroElecCrossSectionDataSet::LoadData(const G4String & argFileName)
{
  CleanUpComponents();

  G4String fullFileName(FullFileName(argFileName));
  std::ifstream in(fullFileName, std::ifstream::binary|std::ifstream::in);

  if (!in.is_open())
    {
      G4String message("Data file \"");
      message+=fullFileName;
      message+="\" not found";
      G4Exception("G4MicroElecCrossSectionDataSet::LoadData","em0003",
                      FatalException,message);
      return false;
    }

  std::vector<G4DataVector *> columns;
  std::vector<G4DataVector *> log_columns;

  std::stringstream *stream(new std::stringstream);
  char c;
  G4bool comment(false);
  G4bool space(true);
  G4bool first(true);

  try
    {
      while (!in.eof())
	{
	  in.get(c);
   
	  switch (c)
	    {
	    case '\r':
	    case '\n':
	      if (!first)
		{
		  unsigned long i(0);
		  G4double value;
      
		  while (!stream->eof())
		    {
		      (*stream) >> value;
       
		      while (i>=columns.size())
                        {
			 columns.push_back(new G4DataVector);
                         log_columns.push_back(new G4DataVector);
                        }
      
		      columns[i]->push_back(value);

// N. A. Karakatsanis
// A condition is applied to check if negative or zero values are present in the dataset.
// If yes, then a near-zero value is applied to allow the computation of the logarithmic value
// If a value is zero, this simplification is acceptable
// If a value is negative, then it is not acceptable and the data of the particular column of
// logarithmic values should not be used by interpolation methods.
//
// Therefore, G4LogLogInterpolation and G4LinLogLogInterpolation should not be used if negative values are present.
// Instead, G4LinInterpolation is safe in every case
// SemiLogInterpolation is safe only if the energy columns are non-negative
// G4LinLogInterpolation is safe only if the cross section data columns are non-negative

                      if (value <=0.) value = 1e-300;
                      log_columns[i]->push_back(std::log10(value));
       
		      i++;
		    }
      
		  delete stream;
		  stream=new std::stringstream;
		}
     
	      first=true;
	      comment=false;
	      space=true;
	      break;

	    case '#':
	      comment=true;
	      break;
     
	    case '\t':
	    case ' ':
	      space = true;
              break;

	    default:
	      if (comment) { break; }
	      if (space && (!first)) { (*stream) << ' '; }

	      first=false;
	      (*stream) << c;
	      space=false;
	    }
	}
    }
  catch(const std::ios::failure &e)
    {
      // some implementations of STL could throw a "failture" exception
      // when read wants read characters after end of file
    }
 
  delete stream;
 
  std::vector<G4DataVector *>::size_type maxI(columns.size());
 
  if (maxI<2)
    {
      G4String message("Data file \"");
      message+=fullFileName;
      message+="\" should have at least two columns";
      G4Exception("G4MicroElecCrossSectionDataSet::LoadData","em0005",
                      FatalException,message);
      return false;
    }
 
  std::vector<G4DataVector*>::size_type i(1);
  while (i<maxI)
    {
      G4DataVector::size_type maxJ(columns[i]->size());

      if (maxJ!=columns[0]->size())
	{
	  G4String message("Data file \"");
	  message+=fullFileName;
	  message+="\" has lines with a different number of columns";
	  G4Exception("G4MicroElecCrossSectionDataSet::LoadData","em0005",
                      FatalException,message);
	  return false;
	}

      G4DataVector::size_type j(0);

      G4DataVector *argEnergies=new G4DataVector;
      G4DataVector *argData=new G4DataVector;
      G4DataVector *argLogEnergies=new G4DataVector;
      G4DataVector *argLogData=new G4DataVector;

      while(j<maxJ)
	{
	  argEnergies->push_back(columns[0]->operator[] (j)*GetUnitEnergies());
	  argData->push_back(columns[i]->operator[] (j)*GetUnitData());
	  argLogEnergies->push_back(log_columns[0]->operator[] (j) + std::log10(GetUnitEnergies()));
	  argLogData->push_back(log_columns[i]->operator[] (j) + std::log10(GetUnitData()));
	  j++;
	}

      AddComponent(new G4EMDataSet(i-1, argEnergies, argData, argLogEnergies, argLogData, GetAlgorithm()->Clone(), GetUnitEnergies(), GetUnitData()));
  
      i++;
    }

  i=maxI;
  while (i>0)
    {
      i--;
      delete columns[i];
      delete log_columns[i];
    }

  return true;
}


G4bool G4MicroElecCrossSectionDataSet::LoadNonLogData(const G4String & argFileName)
{
  CleanUpComponents();

  G4String fullFileName(FullFileName(argFileName));
  std::ifstream in(fullFileName, std::ifstream::binary|std::ifstream::in);

  if (!in.is_open())
    {
      G4String message("Data file \"");
      message+=fullFileName;
      message+="\" not found";
      G4Exception("G4MicroElecCrossSectionDataSet::LoadData","em0003",
                      FatalException,message);
      return false;
    }

  std::vector<G4DataVector *> columns;

  std::stringstream *stream(new std::stringstream);
  char c;
  G4bool comment(false);
  G4bool space(true);
  G4bool first(true);

  try
    {
      while (!in.eof())
	{
	  in.get(c);
   
	  switch (c)
	    {
	    case '\r':
	    case '\n':
	      if (!first)
		{
		  unsigned long i(0);
		  G4double value;
      
		  while (!stream->eof())
		    {
		      (*stream) >> value;
       
		      while (i>=columns.size())
                        {
			 columns.push_back(new G4DataVector);
                        }
      
		      columns[i]->push_back(value);
       
		      i++;
		    }
      
		  delete stream;
		  stream=new std::stringstream;
		}
     
	      first=true;
	      comment=false;
	      space=true;
	      break;
     
	    case '#':
	      comment=true;
	      break;
     
	    case '\t':
	    case ' ':
	      space = true;
              break;

	    default:
	      if (comment) { break; }
	      if (space && (!first)) { (*stream) << ' '; }

	      first=false;
	      (*stream) << c;
	      space=false;
	    }
	}
    }
  catch(const std::ios::failure &e)
    {
      // some implementations of STL could throw a "failture" exception
      // when read wants read characters after end of file
    }
 
  delete stream;
 
  std::vector<G4DataVector *>::size_type maxI(columns.size());
 
  if (maxI<2)
    {
      G4String message("Data file \"");
      message+=fullFileName;
      message+="\" should have at least two columns";
      G4Exception("G4MicroElecCrossSectionDataSet::LoadData","em0005",
                      FatalException,message);
      return false;
    }
 
  std::vector<G4DataVector*>::size_type i(1);
  while (i<maxI)
    {
      G4DataVector::size_type maxJ(columns[i]->size());

      if (maxJ!=columns[0]->size())
	{
	  G4String message("Data file \"");
	  message+=fullFileName;
	  message+="\" has lines with a different number of columns.";
          G4Exception("G4MicroElecCrossSectionDataSet::LoadData","em0005",
                      FatalException,message);
	  return false;
	}

      G4DataVector::size_type j(0);

      G4DataVector *argEnergies=new G4DataVector;
      G4DataVector *argData=new G4DataVector;

      while(j<maxJ)
	{
	  argEnergies->push_back(columns[0]->operator[] (j)*GetUnitEnergies());
	  argData->push_back(columns[i]->operator[] (j)*GetUnitData());
	  j++;
	}

      AddComponent(new G4EMDataSet(i-1, argEnergies, argData, GetAlgorithm()->Clone(), GetUnitEnergies(), GetUnitData()));
  
      i++;
    }

  i=maxI;
  while (i>0)
    {
      i--;
      delete columns[i];
    }

  return true;
}


G4bool G4MicroElecCrossSectionDataSet::SaveData(const G4String & argFileName) const
{
  const size_t n(NumberOfComponents());
 
  if (n==0)
    {
      G4Exception("G4MicroElecCrossSectionDataSet::SaveData","em0005",
                      FatalException,"Expected at least one component");

      return false;
    }
 
  G4String fullFileName(FullFileName(argFileName));
  std::ofstream out(fullFileName);

  if (!out.is_open())
    {
      G4String message("Cannot open \"");
      message+=fullFileName;
      message+="\"";
      G4Exception("G4MicroElecCrossSectionDataSet::SaveData","em0005",
                      FatalException,message);
      return false;
    }
 
  G4DataVector::const_iterator iEnergies(GetComponent(0)->GetEnergies(0).begin());
  G4DataVector::const_iterator iEnergiesEnd(GetComponent(0)->GetEnergies(0).end());
  G4DataVector::const_iterator * iData(new G4DataVector::const_iterator[n]);
 
  size_t k(n);
 
  while (k>0)
    {
      k--;
      iData[k]=GetComponent(k)->GetData(0).begin();
    }
 
  while (iEnergies!=iEnergiesEnd)
    {
      out.precision(10);
      out.width(15);
      out.setf(std::ofstream::left);
      out << ((*iEnergies)/GetUnitEnergies());
  
      k=0;
  
      while (k<n)
	{
	  out << ' ';
	  out.precision(10);
	  out.width(15);
	  out.setf(std::ofstream::left);
	  out << ((*(iData[k]))/GetUnitData());

	  iData[k]++;
	  k++;
	}
  
      out << std::endl;
  
      iEnergies++;
    }
 
  delete[] iData;

  return true;
}


G4String G4MicroElecCrossSectionDataSet::FullFileName(const G4String& argFileName) const
{
  char* path = getenv("G4LEDATA");
  if (!path)
  {
      G4Exception("G4MicroElecCrossSectionDataSet::FullFileName","em0006",
                      FatalException,"G4LEDATA environment variable not set.");
      
      return "";
  }
  
  std::ostringstream fullFileName;
 
  fullFileName << path << "/" << argFileName << ".dat";
                      
  return G4String(fullFileName.str().c_str());
}


G4double G4MicroElecCrossSectionDataSet::FindValue(G4double argEnergy, G4int /* argComponentId */) const
{
  // Returns the sum over the shells corresponding to e
  G4double value = 0.;

  std::vector<G4VEMDataSet *>::const_iterator i(components.begin());
  std::vector<G4VEMDataSet *>::const_iterator end(components.end());

  while (i!=end)
    {
      value+=(*i)->FindValue(argEnergy);
      i++;
    }

  return value;
}


void G4MicroElecCrossSectionDataSet::PrintData(void) const
{
  const size_t n(NumberOfComponents());

  G4cout << "The data set has " << n << " components" << G4endl;
  G4cout << G4endl;
 
  size_t i(0);
 
  while (i<n)
    {
      G4cout << "--- Component " << i << " ---" << G4endl;
      GetComponent(i)->PrintData();
      i++;
    }
}


void G4MicroElecCrossSectionDataSet::SetEnergiesData(G4DataVector* argEnergies, 
					 G4DataVector* argData, 
					 G4int argComponentId)
{
  G4VEMDataSet * component(components[argComponentId]);
 
  if (component)
    {
      component->SetEnergiesData(argEnergies, argData, 0);
      return;
    }

  std::ostringstream message;
  message << "Component " << argComponentId << " not found";
 
  G4Exception("G4MicroElecCrossSectionDataSet::SetEnergiesData","em0005",
                 FatalException,message.str().c_str());
  
}


void G4MicroElecCrossSectionDataSet::SetLogEnergiesData(G4DataVector* argEnergies, 
					 G4DataVector* argData,
                                         G4DataVector* argLogEnergies,
                                         G4DataVector* argLogData, 
					 G4int argComponentId)
{
  G4VEMDataSet * component(components[argComponentId]);
 
  if (component)
    {
      component->SetLogEnergiesData(argEnergies, argData, argLogEnergies, argLogData, 0);
      return;
    }

  std::ostringstream message;
  message << "Component " << argComponentId << " not found";
 
  G4Exception("G4MicroElecCrossSectionDataSet::SetLogEnergiesData","em0005",
                 FatalException,message.str().c_str());

}


void G4MicroElecCrossSectionDataSet::CleanUpComponents()
{
  while (!components.empty())
    {
      if (components.back()) delete components.back();
      components.pop_back();
    }
}


