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
// $Id: G4DNACrossSectionDataSet.cc,v 1.1.2.1 2006/06/29 19:39:02 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// Author: Riccardo Capra <capra@ge.infn.it>
//
// History:
// -----------
// 30 Jun 2005  RC         Created


#include "G4DNACrossSectionDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4EMDataSet.hh"
#include <vector>
#include <sstream>


                                                G4DNACrossSectionDataSet :: G4DNACrossSectionDataSet(G4VDataSetAlgorithm* argAlgorithm, G4double argUnitEnergies, G4double argUnitData)
:
 G4ShellEMDataSet(0, argAlgorithm, argUnitEnergies, argUnitData)
{
}



                                                G4DNACrossSectionDataSet :: ~G4DNACrossSectionDataSet()
{
}





G4bool                                          G4DNACrossSectionDataSet :: LoadData(const G4String & argFileName)
{
 CleanUpComponents();

 G4String fullFileName(FullFileName(argFileName));
 std::ifstream in(fullFileName, std::ifstream::binary|std::ifstream::in);

 if (!in.is_open())
 {
  G4String message("G4DNACrossSectionDataSet::LoadData - data file \"");
  message+=fullFileName;
  message+="\" not found";
  G4Exception(message);
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
        columns.push_back(new G4DataVector);
      
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
     c=' ';
    case ' ':
     if (space)
      break;
    default:
     if (comment)
      break;
     
     if (c==' ')
      space=true;
     else
     {
      if (space && (!first))
       (*stream) << ' ';
      
      first=false;
      (*stream) << c;
      space=false;
     }
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
  G4String message("G4DNACrossSectionDataSet::LoadData - data file \"");
  message+=fullFileName;
  message+="\" should have at least two columns";
  G4Exception(message);
  return false;
 }
 
 std::vector<G4DataVector *>::size_type i(1);
 while (i<maxI)
 {
  G4DataVector::size_type maxJ(columns[i]->size());

  if (maxJ!=columns[0]->size())
  {
   G4String message("G4DNACrossSectionDataSet::LoadData - data file \"");
   message+=fullFileName;
   message+="\" has lines with a different number of columns";
   G4Exception(message);
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



G4bool                                          G4DNACrossSectionDataSet :: SaveData(const G4String & argFileName) const
{
 const size_t n(NumberOfComponents());
 
 if (n==0)
 {
  G4Exception("G4EMDataSet::SaveData - expected at least one component");
  return false;
 }
 
 G4String fullFileName(FullFileName(argFileName));
 std::ofstream out(fullFileName);

 if (!out.is_open())
 {
  G4String message("G4EMDataSet::SaveData - cannot open \"");
  message+=fullFileName;
  message+="\"";
  G4Exception(message);
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





G4String                                        G4DNACrossSectionDataSet :: FullFileName(const G4String & argFileName) const
{
 char* path = getenv("G4LEDATA");
 if (!path)
  G4Exception("G4DNACrossSectionDataSet::FullFileName - G4LEDATA environment variable not set");
  
 std::ostringstream fullFileName;
 
 fullFileName << path << '/' << argFileName << ".dat";
                      
 return G4String(fullFileName.str().c_str());
}
