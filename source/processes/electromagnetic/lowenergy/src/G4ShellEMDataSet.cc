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
// $Id: G4ShellEMDataSet.cc,v 1.14 2006/06/29 19:41:23 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 1 Aug 2001   MGP        Created
// 09.10.01   V.Ivanchenko Add case z=0
//
// -------------------------------------------------------------------

#include "G4ShellEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include <fstream>
#include <sstream>


                                                G4ShellEMDataSet :: G4ShellEMDataSet(G4int argZ, G4VDataSetAlgorithm* argAlgorithm, G4double argUnitEnergies, G4double argUnitData)
:
 z(argZ),
 algorithm(argAlgorithm),
 unitEnergies(argUnitEnergies),
 unitData(argUnitData)
{
 if (algorithm == 0) 
  G4Exception("G4ShellEMDataSet::G4ShellEMDataSet - interpolation == 0");
}



                                                G4ShellEMDataSet :: ~G4ShellEMDataSet()
{
 CleanUpComponents();
 
 if (algorithm)
  delete algorithm;
}





G4double                                        G4ShellEMDataSet :: FindValue(G4double argEnergy, G4int /* argComponentId */) const
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





void                                            G4ShellEMDataSet :: PrintData(void) const
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





void                                            G4ShellEMDataSet :: SetEnergiesData(G4DataVector * argEnergies, G4DataVector * argData, G4int argComponentId)
{
 G4VEMDataSet * component(components[argComponentId]);
 
 if (component)
 {
  component->SetEnergiesData(argEnergies, argData, 0);
  return;
 }

 std::ostringstream message;
 message << "G4ShellEMDataSet::SetEnergiesData - component " << argComponentId << " not found";
 
 G4Exception(message.str().c_str());
}





G4bool                                          G4ShellEMDataSet :: LoadData(const G4String & argFileName)
{
 CleanUpComponents();

 G4String fullFileName(FullFileName(argFileName));
 std::ifstream in(fullFileName);

 if (!in.is_open())
 {
  G4String message("G4ShellEMDataSet::LoadData - data file \"");
  message+=fullFileName;
  message+="\" not found";
  G4Exception(message);
 }

 G4DataVector * argEnergies(0);
 G4DataVector * argData(0);

 G4double a;
 G4int shellIndex(0);
 bool energyColumn(true);

 do
 {
  in >> a;
  
  if (a == -1)
  {
   if (energyColumn && argEnergies!=0)
   {
    AddComponent(new G4EMDataSet(shellIndex, argEnergies, argData, algorithm->Clone(), unitEnergies, unitData));
    argEnergies=0;
    argData=0;
   }
   
   energyColumn=(!energyColumn);
  }
  else if (a!=-2)
  {
   if (argEnergies==0)
   {
    argEnergies=new G4DataVector;
    argData=new G4DataVector;
   }
  
   if (energyColumn)
    argEnergies->push_back(a*unitEnergies);
   else
    argData->push_back(a*unitData);

   energyColumn=(!energyColumn);
  }
 }
 while (a != -2);

 return true;
}



G4bool                                          G4ShellEMDataSet :: SaveData(const G4String & argFileName) const
{
 G4String fullFileName(FullFileName(argFileName));
 std::ofstream out(fullFileName);

 if (!out.is_open())
 {
  G4String message("G4EMDataSet::SaveData - cannot open \"");
  message+=fullFileName;
  message+="\"";
  G4Exception(message);
 }
 
 const size_t n(NumberOfComponents());
 size_t k(0);
 
 while (k<n)
 {
  const G4VEMDataSet * component=GetComponent(k);
  
  if (component)
  {
   const G4DataVector & energies(component->GetEnergies(0));
   const G4DataVector & data(component->GetData(0));
 
   G4DataVector::const_iterator i(energies.begin());
   G4DataVector::const_iterator endI(energies.end());
   G4DataVector::const_iterator j(data.begin());
  
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





void                                            G4ShellEMDataSet :: CleanUpComponents(void)
{
 while (!components.empty())
 {
  if (components.back())
   delete components.back();

  components.pop_back();
 }
}





G4String                                        G4ShellEMDataSet :: FullFileName(const G4String & argFileName) const
{
 char* path = getenv("G4LEDATA");
 if (!path)
  G4Exception("G4ShellEMDataSet::FullFileName - G4LEDATA environment variable not set");
  
 std::ostringstream fullFileName;
 
 fullFileName << path << '/' << argFileName << z << ".dat";
                      
 return G4String(fullFileName.str().c_str());
}
