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
// $Id: G4ShellVacancyTest.cc,v 1.1 2001-09-21 10:23:23 elena Exp $
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
#include "g4std/fstream"
#include "g4std/iomanip"
#include "G4CompositeEMDataSet.hh"
#include "G4ShellEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4VEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4ShellVacancy.hh"
#include "G4Material.hh"

int main()
{ 
  //--------- Materials definition ---------

  G4Material* Be = new G4Material("Beryllium",    4.,  9.01*g/mole, 1.848*g/cm3);
  G4Material* Graphite = new G4Material("Graphite",6., 12.00*g/mole, 2.265*g/cm3 );
  G4Material* Al  = new G4Material("Aluminium", 13., 26.98*g/mole, 2.7 *g/cm3);
  //G4Material* Si  = new G4Material("Silicon",   14., 28.055*g/mole, 2.33*g/cm3);
  //G4Material* LAr = new G4Material("LArgon",   18., 39.95*g/mole, 1.393*g/cm3);
  //G4Material* Fe  = new G4Material("Iron",      26., 55.85*g/mole, 7.87*g/cm3);
  //G4Material* Cu  = new G4Material("Copper",    29., 63.55*g/mole, 8.96*g/cm3);
  //G4Material*  W  = new G4Material("Tungsten", 74., 183.85*g/mole, 19.30*g/cm3);
  //G4Material* Pb  = new G4Material("Lead",      82., 207.19*g/mole, 11.35*g/cm3);
  //G4Material*  U  = new G4Material("Uranium", 92., 238.03*g/mole, 18.95*g/cm3);

  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I  = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);

  //G4Material*  maO = new G4Material("Oxygen", 8., 16.00*g/mole, 1.1*g/cm3);

  G4Material* water = new G4Material ("Water" , 1.*g/cm3, 2);
  water->AddElement(H,2);
  water->AddElement(O,1);

  G4Material* ethane = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  ethane->AddElement(H,6);
  ethane->AddElement(C,2);
  
  G4Material* csi = new G4Material ("CsI" , 4.53*g/cm3, 2);
  csi->AddElement(Cs,1);
  csi->AddElement(I,1);

  G4ShellVacancy* manager = new G4ShellVacancy;  
// Setup

  G4cout << "G4EMDataSet test: dump LLNL data sets" << G4endl;

  G4cout.setf( ios::scientific, ios::floatfield );

  G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();

  G4cout << "Interpolation created" << G4endl; 
    G4String  fileName = "brem/br-cs-";

    G4VEMDataSet* dataSet1;

      dataSet1 = new G4CompositeEMDataSet(fileName,interpolation,MeV,10000);
    
      fileName = "comp/ce-cs-";
	
      G4VEMDataSet* dataSet2;
      dataSet2 = new G4CompositeEMDataSet(fileName,interpolation,MeV,10000);
     
      fileName = "comp/ce-sf-";
      G4VEMDataSet* dataSet3;

      dataSet3 = new G4CompositeEMDataSet(fileName,interpolation,MeV,10000);


      fileName = "pair/pp-cs-";
      G4VEMDataSet* dataSet4;

      dataSet4 = new G4CompositeEMDataSet(fileName,interpolation,MeV,10000);
      static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

      fileName = "phot/pe-cs-";
      G4VEMDataSet* dataSet5;

      dataSet5 = new G4CompositeEMDataSet(fileName,interpolation,MeV,10000);

      fileName = "phot/pe-ss-cs-";
      G4VEMDataSet* dataSet6;

      dataSet6 = new G4CompositeEMDataSet(fileName,interpolation,MeV,10000);

      G4int nMaterials = theMaterialTable->length();
      
      G4cout << "Available materials are: " << G4endl;
      for (G4int mat = 0; mat < nMaterials; mat++)
    {
      G4cout << mat << ") "
	     << (*theMaterialTable)(mat)->GetName()
	     << G4endl;
    }
      G4int materialId=-1;
      G4cout << "Which material? " << G4endl;
      G4cin >> materialId;

      G4Material* material = (*theMaterialTable)(materialId) ;
      
      G4cout << "The selected material is: "
	 << material->GetName()
	     << G4endl;

 

  G4cout << "Enter incidentEnergy" << G4endl;
  G4double e;
  G4cin >> e;

  G4cout << "Enter energyLoss" << G4endl;
  G4double eLoss;
  G4cin >> eLoss;

  manager->AddXsiTable(dataSet1);
  manager->AddXsiTable(dataSet2);
  manager->AddXsiTable(dataSet3);
  manager->AddXsiTable(dataSet4);
  manager->AddXsiTable(dataSet5);
  manager->AddXsiTable(dataSet6);

  G4std::vector<G4int> vector = manager->GenerateNumberOfIonisations(material,
								    e,eLoss);				    
  size_t vectorSize = vector.size();
  for (size_t p=0; p<vectorSize; p++)
    {
      G4int n = vector[p];
      G4cout<<"The number of transition for the "<<p<<"th element in the material "<<
	material->GetName()<<" is : "<<n<<G4endl;
    }

  delete manager;

  cout << "END OF THE MAIN PROGRAM" << G4endl;
}












