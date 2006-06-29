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
// $Id: G4EMDataSetTest.cc,v 1.3 2006-06-29 19:43:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4EMDataSetTest
//
//      Author:        Maria Grazia Pia
// 
//      Creation date: 1 August 2001
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include "G4CompositeEMDataSet.hh"
#include "G4ShellEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4VEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"

int main()
{
  // Setup

  G4cout << "G4EMDataSet test: dump LLNL data sets" << G4endl;

  G4cout << "Enter Z " << G4endl;
  G4int Z;
  G4cin >> Z;

  G4cout << "Enter file name " << G4endl;
  G4String fileName;
  //  G4cin >> fileName;
  fileName = "brem/br-cs-";

  G4cout.setf( ios::scientific, ios::floatfield );

  G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();

  G4cout << "Interpolation created" << G4endl; 

  G4cout << "Atom (1) or Shell (2) or Composite (3) data set?" << G4endl;
  G4int type;
  G4cin >> type;

  G4VEMDataSet* dataSet;

  if (type == 1)  
    {
      fileName = "brem/br-cs-";
      dataSet = new G4EMDataSet(Z,fileName,interpolation);
    }
  else if (type == 2) 
    {
      fileName = "phot/pe-ss-cs-";
      dataSet = new G4ShellEMDataSet(Z,fileName,interpolation);
    }
  else
    {
      fileName = "brem/br-cs-";
      dataSet = new G4CompositeEMDataSet(fileName,interpolation);
    }

  dataSet->PrintData();

  G4cout << "Enter energy" << G4endl;
  G4double e;
  G4cin >> e;

  G4double sigma = dataSet->FindValue(e) / barn;
  G4cout << "Value = " << sigma << G4endl;

  G4cout << "Enter shell index" << G4endl;
  G4int id;
  G4cin >> id;

  const G4VEMDataSet* component = dataSet->GetComponent(id);
  if (component != 0)
    {
      sigma = component->FindValue(e) / barn;
      G4cout << "Value for id = " << id << ": " << sigma << G4endl;
      G4DataVector data = component->GetData(id);
      G4int sizeD = data.size();
      G4cout << "The energy vector of the component has " << sizeD << " elements" << G4endl;
    }

  G4DataVector energies = dataSet->GetEnergies(id);
  G4int sizeE = energies.size();
  G4cout << "The energy vector of component " << id 
	 << " has " << sizeE << " elements" << G4endl;

  delete dataSet;

  cout << "END OF THE MAIN PROGRAM" << G4endl;
}












