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
// $Id: G4CrossSectionHandlerTest.cc,v 1.1 2001-08-20 17:28:50 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4DataHandlerTest
//
//      Author:        Maria Grazia Pia
// 
//      Creation date: 3 August 2001
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4CrossSectionHandler.hh"
#include "G4VEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"

HepTupleManager* hbookManager;

int main()
{
  G4cout.setf( ios::scientific, ios::floatfield );


 // ---- HBOOK initialization


  hbookManager = new HBookFile("datahandler.hbook", 58);
  assert (hbookManager != 0);
  
  // ---- primary ntuple ------
  HepTuple* ntuple = hbookManager->ntuple("CrossSectionNtuple");
  assert (ntuple != 0);
  
  // Create some materials 
  G4Material* Be = new G4Material("Beryllium",    4.,  9.01*g/mole, 1.848*g/cm3);
  G4Material* Graphite = new G4Material("Graphite",6., 12.00*g/mole, 2.265*g/cm3 );
  G4Material* Al  = new G4Material("Aluminium", 13., 26.98*g/mole, 2.7 *g/cm3);
  G4Material* Si  = new G4Material("Silicon",   14., 28.055*g/mole, 2.33*g/cm3);
  G4Material* LAr = new G4Material("LArgon",   18., 39.95*g/mole, 1.393*g/cm3);
  G4Material* Fe  = new G4Material("Iron",      26., 55.85*g/mole, 7.87*g/cm3);
  G4Material* Cu  = new G4Material("Copper",    29., 63.55*g/mole, 8.96*g/cm3);
  G4Material*  W  = new G4Material("Tungsten", 74., 183.85*g/mole, 19.30*g/cm3);
  G4Material* Pb  = new G4Material("Lead",      82., 207.19*g/mole, 11.35*g/cm3);
  G4Material*  U  = new G4Material("Uranium", 92., 238.03*g/mole, 18.95*g/cm3);

  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I  = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);

  G4Material*  maO = new G4Material("Oxygen", 8., 16.00*g/mole, 1.1*g/cm3);

  G4Material* water = new G4Material ("Water" , 1.*g/cm3, 2);
  water->AddElement(H,2);
  water->AddElement(O,1);

  G4Material* ethane = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  ethane->AddElement(H,6);
  ethane->AddElement(C,2);
  
  G4Material* csi = new G4Material ("CsI" , 4.53*g/cm3, 2);
  csi->AddElement(Cs,1);
  csi->AddElement(I,1);

  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4int nMaterials = theMaterialTable->length();
  G4cout << "The MaterialTable contains "
	 << nMaterials
	 << " materials "
	 << G4endl;

  G4cout << "G4DataSetManager test: dump LLNL data sets" << G4endl;

  //  G4cout << "Enter file name " << G4endl;

  G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();
  G4cout << "Interpolation created" << G4endl;

  G4CrossSectionHandler*  manager = new G4CrossSectionHandler(interpolation,10*eV);
  G4cout << "CrossSectionHandler created" << G4endl;
 
  G4cout << " Load Atom (1) or Shell (2) data?" << G4endl;
  G4int type;
  G4cin >> type;
  G4String fileName;
  if (type == 1)
    {
      // G4cin >> fileName;
      fileName = "brem/br-cs-";
      G4String name(fileName);
      manager->LoadData(fileName);
      G4cout << "Data Loaded" << G4endl;
    }
  if (type == 2)
    {
      // G4cin >> fileName;
      fileName = "phot/pe-ss-cs-";
      G4String name(fileName);
      manager->LoadShellData(fileName);
      G4cout << "Data Loaded" << G4endl;
    }
  
  G4cout << "Data loaded; print (1) or continue (2)? " << G4endl;
  G4int printData;
  G4cin >> printData;
  if (printData == 1)  manager->PrintData();

  G4cout << "Enter Z " << G4endl;
  G4int Z;
  G4cin >> Z;

  G4cout << "Enter energy" << G4endl;
  G4double e;
  G4cin >> e;
  G4double sigma = manager->FindValue(Z,e) / barn;
  G4cout << "Value retrieved by manager = " << sigma << G4endl;

  G4cout << "Available materials are: " << G4endl;
  for (G4int mat = 0; mat < nMaterials; mat++)
    {
      G4cout << mat << ") "
	     << (*theMaterialTable)(mat)->GetName()
	     << G4endl;
    }

  G4int materialId;
  G4cout << "Which material? " << G4endl;
  G4cin >> materialId;
  G4Material* material = (*theMaterialTable)(materialId) ;

  G4double materialSigma = manager->ValueForMaterial(material,e);
  G4cout << "Material value calculated by manager = " << materialSigma << G4endl;

  G4int selectedZ = manager->SelectRandomAtom(material,e);
  G4cout << "Select random atom: Z = " << selectedZ << G4endl;

  G4int selectedShell = manager->SelectRandomShell(Z,e);
  G4cout << "Select random shell: shell = " << selectedShell << G4endl;

  G4cout << "Now clearing" << G4endl;
  manager->Clear();

  G4cout << "Now re-loading" << G4endl;
  if (type == 1) manager->LoadData(fileName);
  if (type == 2) manager->LoadShellData(fileName);
  G4cout << "End of re-loading" << G4endl;

  G4VEMDataSet* meanFreePathTable = manager->BuildMeanFreePathForMaterials();
  size_t nRows = meanFreePathTable->NumberOfComponents();
  meanFreePathTable->PrintData();
  G4cout << "MeanFreePathTable has " << nRows << " components" << G4endl;

  const G4VEMDataSet* materialSet = meanFreePathTable->GetComponent(materialId);

  G4int i;
  for (i=2; i<10; i++)
    {
      G4double e = (10. * i) * eV;
      G4double value = manager->FindValue(Z,e) / barn;
      G4double valueM = manager->ValueForMaterial(material,e);
      G4double valueT = materialSet->FindValue(e);
      if (valueT > 0.) valueT = 1. / valueT;
      ntuple->column("e",e);
      ntuple->column("sigmac",valueM);
      ntuple->column("sigmat",valueT);
      ntuple->dumpData();
      
      //     G4cout << "e = " << e/MeV << " --- Cross section = " << value << " b" << G4endl;
    }
  for (i=1; i<10; i++)
    {
      G4double e = (100. * i) * eV;
      G4double value = manager->FindValue(Z,e) / barn;
      //    G4cout << "e = " << e/MeV << " --- Cross section = " << value << " b" << G4endl;
      G4double valueM = manager->ValueForMaterial(material,e);
      G4double valueT = materialSet->FindValue(e);
      if (valueT > 0.) valueT = 1. / valueT;
      ntuple->column("e",e);
      ntuple->column("sigmac",valueM);
      ntuple->column("sigmat",valueT);
      ntuple->dumpData();  
    }
  for (i=1; i<10; i++)
    {
      G4double e = (1. * i) * keV;
      G4double value = manager->FindValue(Z,e) / barn;
      G4double valueM = manager->ValueForMaterial(material,e);
      G4double valueT = materialSet->FindValue(e);
      if (valueT > 0.) valueT = 1. / valueT;
      ntuple->column("e",e);
      ntuple->column("sigmac",valueM);
      ntuple->column("sigmat",valueT);
      ntuple->dumpData();  
      //     G4cout << "e = " << e/MeV << " --- Cross section = " << value << " b" << G4endl;
    }
  for (i=1; i<10; i++)
    {
      G4double e = (10. * i) * keV;
      G4double value = manager->FindValue(Z,e) / barn;
      G4double valueM = manager->ValueForMaterial(material,e);
      G4double valueT = materialSet->FindValue(e);
      if (valueT > 0.) valueT = 1. / valueT;
      ntuple->column("e",e);
      ntuple->column("sigmac",valueM);
      ntuple->column("sigmat",valueT);
      ntuple->dumpData();
      
      //      G4cout << "e = " << e/MeV << " --- Cross section = " << value << " b" << G4endl;
    }
  for (i=1; i<10; i++)
    {
      G4double e = (100. * i) * keV;
      G4double value = manager->FindValue(Z,e) / barn;
      G4double valueM = manager->ValueForMaterial(material,e);
      G4double valueT = materialSet->FindValue(e);
      if (valueT > 0.) valueT = 1. / valueT;
      ntuple->column("e",e);
      ntuple->column("sigmac",valueM);
      ntuple->column("sigmat",valueT);
      ntuple->dumpData();
      
      //      G4cout << "e = " << e/MeV << " --- Cross section = " << value << " b" << G4endl;
    }
  for (i=1; i<10; i++)
    {
      G4double e = (1. * i) * MeV;
      G4double value = manager->FindValue(Z,e) / barn;
      G4double valueM = manager->ValueForMaterial(material,e);
      G4double valueT = materialSet->FindValue(e);
      if (valueT > 0.) valueT = 1. / valueT;
      ntuple->column("e",e);
      ntuple->column("sigmac",valueM);
      ntuple->column("sigmat",valueT);
      ntuple->dumpData();
      
      //      G4cout << "e = " << e/MeV << " --- Cross section = " << value << " b" << G4endl;
    }
  for (i=1; i<10; i++)
    {
      G4double e = (10. * i) * MeV;
      G4double value = manager->FindValue(Z,e) / barn;
      G4double valueM = manager->ValueForMaterial(material,e);
      G4double valueT = materialSet->FindValue(e);
      if (valueT > 0.) valueT = 1. / valueT;
      ntuple->column("e",e);
      ntuple->column("sigmac",valueM);
      ntuple->column("sigmat",valueT);
      ntuple->dumpData();
      //      G4cout << "e = " << e/MeV << " --- Cross section = " << value << " b" << G4endl;
    }
  for (i=1; i<10; i++)
    {
      G4double e = (100. * i) * MeV;
      G4double value = manager->FindValue(Z,e) / barn;
      G4double valueM = manager->ValueForMaterial(material,e);
      G4double valueT = materialSet->FindValue(e);
      if (valueT > 0.) valueT = 1. / valueT;
      ntuple->column("e",e);
      ntuple->column("sigmac",valueM);
      ntuple->column("sigmat",valueT);
      ntuple->dumpData();
      
      //      G4cout << "e = " << e/MeV << " --- Cross section = " << value << " b" << G4endl;
    }
  for (i=1; i<10; i++)
    {
      G4double e = (1. * i) * GeV;
      G4double value = manager->FindValue(Z,e) / barn;
      G4double valueM = manager->ValueForMaterial(material,e);
      G4double valueT = materialSet->FindValue(e);
      if (valueT > 0.) valueT = 1. / valueT;
      ntuple->column("e",e);
      ntuple->column("sigmac",valueM);
      ntuple->column("sigmat",valueT);
      ntuple->dumpData();
      
      //      G4cout << "e = " << e/MeV << " --- Cross section = " << value << " b" << G4endl;
    }
  for (i=1; i<10; i++)
    {
      G4double e = (10. * i) * GeV;
      G4double value = manager->FindValue(Z,e) / barn;
      G4double valueM = manager->ValueForMaterial(material,e);
      G4double valueT = materialSet->FindValue(e);
      if (valueT > 0.) valueT = 1. / valueT;
      ntuple->column("e",e);
      ntuple->column("sigmac",valueM);
      ntuple->column("sigmat",valueT);
      ntuple->dumpData();    
      //      G4cout << "e = " << e/MeV << " --- Cross section = " << value << " b" << G4endl;
    }
  
  G4double eFin = 100 * GeV;
  G4double value = manager->FindValue(Z,eFin) / barn;
  G4double valueM = manager->ValueForMaterial(material,eFin);
  G4double valueT = materialSet->FindValue(eFin);
  if (valueT > 0.) valueT = 1. / valueT;
  ntuple->column("e",eFin);
  ntuple->column("sigmac",valueM);
  ntuple->column("sigmat",valueT);
  ntuple->dumpData();

  //      G4cout << "e = " << e/MeV << " --- Cross section = " << value << " b" << G4endl;

  hbookManager->write();
  delete hbookManager;
  
  delete meanFreePathTable;    
  delete manager;

  cout << "END OF THE MAIN PROGRAM" << G4endl;
}








