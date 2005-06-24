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
//
// $Id: G4EMDataSet.cc,v 1.10 2005-06-24 09:55:05 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "G4EMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include <fstream>
#include <sstream>

                                                 G4EMDataSet :: G4EMDataSet(G4int argZ, const G4VDataSetAlgorithm * argAlgorithm, G4double argUnitEnergies, G4double argUnitData)
:
 z(argZ),
 energies(0),
 data(0),
 algorithm(argAlgorithm),
 unitEnergies(argUnitEnergies),
 unitData(argUnitData)
{
 if (algorithm == 0) 
  G4Exception("G4EMDataSet::G4EMDataSet - interpolation == 0");
}



                                                 G4EMDataSet :: G4EMDataSet(G4int argZ, G4DataVector * argEnergies, G4DataVector * argData, const G4VDataSetAlgorithm * argAlgorithm, G4double argUnitEnergies, G4double argUnitData)
:
 z(argZ),
 energies(argEnergies),
 data(argData),
 algorithm(argAlgorithm),
 unitEnergies(argUnitEnergies),
 unitData(argUnitData)
{
 if (algorithm==0) 
  G4Exception("G4EMDataSet::G4EMDataSet - interpolation == 0");

 if ((energies==0) ^ (data==0))
  G4Exception("G4EMDataSet::G4EMDataSet - different size for energies and data (zero case)");

 if (energies==0)
  return;
  
 if (energies->size()!=data->size())
  G4Exception("G4EMDataSet::G4EMDataSet - different size for energies and data");
}



                                                 G4EMDataSet :: ~G4EMDataSet()
{ 
 delete algorithm;
 
 if (energies)
  delete energies;

 if (data)
  delete data;
}





G4double                                        G4EMDataSet :: FindValue(G4double argEnergy, G4int /* argComponentId */) const
{
 if (!energies)
  G4Exception("G4EMDataSet::FindValue - energies == 0");
  
 if (energies->empty())
  return 0;

 if (argEnergy <= (*energies)[0])
  return (*data)[0];

 size_t i(energies->size()-1);

 if (argEnergy >= (*energies)[i])
  return (*data)[i];

 return algorithm->Calculate(argEnergy, FindLowerBound(argEnergy), *energies, *data);
}





void                                            G4EMDataSet :: PrintData(void) const
{
 if (!energies)
 {
  G4cout << "Data not available." << G4endl;
  return;
 }
  
 size_t size = energies->size();
 
 for (size_t i(0); i<size; i++)
  G4cout << "Point: " << ((*energies)[i]/unitEnergies)
	 << " - Data value : " << ((*data)[i]/unitData) << G4endl; 
}





void                                            G4EMDataSet :: SetEnergiesData(G4DataVector * argEnergies, G4DataVector * argData, G4int /* argComponentId */)
{
 if (energies)
  delete energies;
 energies=argEnergies;

 if (data)
  delete data;
 data=argData;
 
 if ((energies==0) ^ (data==0))
  G4Exception("G4EMDataSet::SetEnergiesData - different size for energies and data (zero case)");

 if (energies==0)
  return;
  
 if (energies->size()!=data->size())
  G4Exception("G4EMDataSet::SetEnergiesData - different size for energies and data");
}





G4bool                                          G4EMDataSet :: LoadData(const G4String & argFileName)
{
 // The file is organized into two columns:
 // 1st column is the energy
 // 2nd column is the corresponding value
 // The file terminates with the pattern: -1   -1
 //                                       -2   -2
 
 G4String fullFileName(FullFileName(argFileName));
 std::ifstream in(fullFileName);

 if (!in.is_open())
 {
  G4String message("G4EMDataSet::LoadData - data file \"");
  message+=fullFileName;
  message+="\" not found";
  G4Exception(message);
 }

 G4DataVector * argEnergies=new G4DataVector;
 G4DataVector * argData=new G4DataVector;

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
 
 return true;
}



G4bool                                          G4EMDataSet :: SaveData(const G4String & argFileName) const
{
 // The file is organized into two columns:
 // 1st column is the energy
 // 2nd column is the corresponding value
 // The file terminates with the pattern: -1   -1
 //                                       -2   -2
 
 G4String fullFileName(FullFileName(argFileName));
 std::ofstream out(fullFileName);

 if (!out.is_open())
 {
  G4String message("G4EMDataSet::SaveData - cannot open \"");
  message+=fullFileName;
  message+="\"";
  G4Exception(message);
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
   out << (*i) << (*j) << std::endl;
   i++;
   j++;
  }
 }
 
 out << -1. << -1. << std::endl << -2. << -2. << std::endl;

 return true;
}





size_t                                          G4EMDataSet :: FindLowerBound(G4double argEnergy) const
{
 size_t lowerBound(0);
 size_t upperBound(energies->size() - 1);
  
 while (lowerBound <= upperBound) 
 {
  size_t midBin((lowerBound + upperBound)/2);

  if (argEnergy < (*energies)[midBin])
   upperBound = midBin-1;
  else
   lowerBound = midBin+1;
 }
  
 return upperBound;
}





G4String                                        G4EMDataSet :: FullFileName(const G4String & argFileName) const
{
 char* path = getenv("G4LEDATA");
 if (!path)
  G4Exception("G4EMDataSet::FullFileName - G4LEDATA environment variable not set");
  
 std::ostringstream fullFileName;
 
 fullFileName << path << '/' << argFileName << z << ".dat";
                      
 return G4String(fullFileName.str().c_str());
}
