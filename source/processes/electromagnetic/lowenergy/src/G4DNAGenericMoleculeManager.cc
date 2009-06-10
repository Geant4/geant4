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
// ----------------- History -----------------
//
//     2 Apr 2009 - A.Mantero - created



#include "G4DNAGenericMoleculeManager.hh"
#include "G4DynamicParticle.hh"

#include "G4Molecule.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAGenericMoleculeManager * G4DNAGenericMoleculeManager :: Instance(void)
{
 if (!theInstance)
  theInstance=new G4DNAGenericMoleculeManager;
 
 return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Molecule* G4DNAGenericMoleculeManager :: GetMolecule(const G4String & name)
{
 MoleculeMap::const_iterator i(map.find(name));
 
 if (i==map.end())
  return 0;
  
 return i->second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAGenericMoleculeManager :: G4DNAGenericMoleculeManager()
{
  
  ionManager = G4DNAGenericIonsManager::Instance();

  // molecules construction
  
  G4Molecule* oxonium; // H3O -- it will become H3O+
  G4Molecule* hydroxyl; // OH -- it will produce OH- too
  G4Molecule* molHydrogen; // H2
  //G4Molecule* hydroxide; // OH-
  G4Molecule* hydroPeroxide; // H2O2
  G4Molecule* water; // H2O -- it will become also H2O+
  G4Molecule* solvElectron; // e_aq it is treated like a molecule: every time you need it, 
                            // you must create it and then add an electron 
  
  // Fisrt I need to have Molecules definition
  // then I can add decay channels, with a proper function.
  // let-s remember how to define a molecule:
  //    G4Molecule(G4String aName, G4double mass, G4String type, G4double lifetime, G4int electronsNumber, G4int electronicLevels, G4double diffCoeff, G4double jumpLenght, G4double diffTSm, G4int atomsNumber, G4double radius);
  

  G4double mass = 19.02*g/Avogadro; 
  
  oxonium = new G4Molecule("H3O", mass, "molecule", 10e-12 * s, 11, 5, 9.0 *10e-9*(m*m/s), 0.23 * nm, 1 * 10e-12 * s, 4, 0.961*2 * angstrom);
  
  oxonium->SetLevelOccupation(0);
  oxonium->SetLevelOccupation(1);
  oxonium->SetLevelOccupation(2,4);
  oxonium->SetLevelOccupation(3);
  oxonium->SetLevelOccupation(4,1); 

  mass = 0; // molecules born neutral: this time we have to do a trick and add an electron when we create the dynamic object.
  
  solvElectron = new G4Molecule("e_aq", mass, "molecule", 1e-12 * s, 0, 1, 4.5 *10e-9*(m*m/s), 0.16 * nm, 1 * 10e-12 * s, 1, 0.23* nm); // radius from K.D. Jordan, Sciencem vol. 306 p.618
  
  solvElectron->SetLevelOccupation(0,0);

  mass = 17.00734*g/Avogadro; 
  
  hydroxyl = new G4Molecule("OH", mass, "molecule", 10e-12 * s, 9, 5, 5.0 *10e-9*(m*m/s), 0.17 * nm, 1 * 10e-12 * s, 2, 0.958 * angstrom);
  
  hydroxyl->SetLevelOccupation(0);
  hydroxyl->SetLevelOccupation(1);
  hydroxyl->SetLevelOccupation(2);
  hydroxyl->SetLevelOccupation(3,3);
  
  
  mass = 2.01588*g/Avogadro;
  
  molHydrogen = new G4Molecule("H2", mass, "molecule", 10e-12 * s, 2, 2, 5.0 *10e-9*(m*m/s), 0.17 * nm, 1 * 10e-12 * s, 2, 0.958 * angstrom);
  
  molHydrogen->SetLevelOccupation(0);
  
  
  mass = 34.01468*g/Avogadro; 
  
  hydroPeroxide = new G4Molecule("H2O2", mass, "molecule", 10e-12 * s, 18, 10, 1.4 *10e-9*(m*m/s), 0.09 * nm, 1 * 10e-12 * s, 2, 3 * angstrom); // some invented values
  
  hydroPeroxide->SetLevelOccupation(0);
  hydroPeroxide->SetLevelOccupation(1);
  hydroPeroxide->SetLevelOccupation(2);
  hydroPeroxide->SetLevelOccupation(3);
  hydroPeroxide->SetLevelOccupation(4);
  hydroPeroxide->SetLevelOccupation(5);
  hydroPeroxide->SetLevelOccupation(6);
  hydroPeroxide->SetLevelOccupation(7);
  
  
  mass = 19.02*g/Avogadro; // this is what we have now
  
  water = new G4Molecule("H2O", mass, "molecule", 10e-12 * s, 10, 8, 0, 0, 0, 3, 2.75 * angstrom);
  
  water->SetLevelOccupation(0);
  water->SetLevelOccupation(1);
  water->SetLevelOccupation(2);
  water->SetLevelOccupation(3);
  water->SetLevelOccupation(4); 
  
  
  map["e_aq"]  =solvElectron;
  map["H3O" ]  =oxonium;
  map["OH"  ]  =hydroxyl;
  map["H2"  ]  =molHydrogen;
  map["H2O2"]  =hydroPeroxide;
  map["H2O" ]  =water; 

  
  CreateDecayChannels();
  
}


void G4DNAGenericMoleculeManager::CreateDecayChannels() {
  //
  // ******************************************
  // *              Building WATER            *
  // ******************************************
  //
  
  G4Molecule* water = GetMolecule("H2O");
  
  G4MolecularDecayChannel* decCh1 = new G4MolecularDecayChannel();
  G4MolecularDecayChannel* decCh2 = new G4MolecularDecayChannel();

  G4DynamicParticle product = G4DynamicParticle();
  G4DynamicParticle product2 = G4DynamicParticle();
  G4DynamicParticle product3 = G4DynamicParticle();

  product.SetDefinition( (G4ParticleDefinition*)GetMolecule("H2O") );

  decCh1->AddProduct(product);//it should work.
  decCh1->SetEnergy(5*eV);
  decCh1->SetProbability(0.35);

  product.SetDefinition( (G4ParticleDefinition*)GetMolecule("OH") );
  product2.SetDefinition(ionManager->GetIon("hydrogen"));

  decCh2->AddProduct(product);//it should work.
  decCh2->AddProduct(product2);//it should work.
  decCh2->SetProbability(0.65);
  
  water->AddExcitedState("A^1B_1");
  
  water->AddDecayChannel("A^1B_1",*decCh1);
  water->AddDecayChannel("A^1B_1",*decCh2);
  
  G4ElectronOccupancy* occ1 = new G4ElectronOccupancy();
  G4ElectronOccupancy* occ2 = new G4ElectronOccupancy();
  G4ElectronOccupancy* occ3 = new G4ElectronOccupancy();
  
  *occ1 = *occ2 = *occ3 = *(water->GetGroundState());   
  
  // this will be the methoid probably used by processes to identify the vacancy.
  
  occ1->RemoveElectron(4); // this is the transition form ground state to 
  occ1->AddElectron(5,1);  // the first unoccupied orbital: A^1B_1        
  
  water->AddeConfToExcitedState("A^1B_1", *occ1);
  
  delete decCh1;
  delete decCh2;
  
  decCh1 = new G4MolecularDecayChannel();
  decCh2 = new G4MolecularDecayChannel();
  G4MolecularDecayChannel* decCh3 = new G4MolecularDecayChannel();

  product.SetDefinition( (G4ParticleDefinition*)GetMolecule("H2O") );

  
  decCh1->AddProduct(product);//it should work.
  decCh1->SetEnergy(5*eV);
  decCh1->SetProbability(0.52);

  product.SetDefinition( (G4ParticleDefinition*)GetMolecule("H2") );
  product.SetDefinition( (G4ParticleDefinition*)GetMolecule("OH") );
  
  decCh2->AddProduct(product);//it should work.
  decCh2->AddProduct(product2);//it should work.
  decCh2->AddProduct(product2);//it should work.
  decCh2->SetProbability(0.145);
  

  product.SetDefinition( (G4ParticleDefinition*)GetMolecule("H3O") );
  product.RemoveElectron(5); // I remove the e- of the 5th level, that is the LUMO of H3O+

  product2.SetDefinition( (G4ParticleDefinition*)GetMolecule("OH") );
  product3.SetDefinition( (G4ParticleDefinition*)GetMolecule("e_aq") );
  product3.AddElectron(0); // Creating the solvated electron.
  
  decCh3->AddProduct(product);
  decCh3->AddProduct(product2);//it should work.
  decCh3->AddProduct(product3);//it should work.
  decCh3->SetProbability(0.335);
  
  occ2->RemoveElectron(3); // this is the transition form ground state to 
  occ2->AddElectron(5,1);  // the first unoccupied orbital: B^1A_1        
  
  water->AddeConfToExcitedState("B^1A_1", *occ2);
  water->AddDecayChannel("A^1B_1",*decCh1);
  water->AddDecayChannel("B^1A_1",*decCh2);
  water->AddDecayChannel("B^1A_1",*decCh3);
  
  occ3->RemoveElectron(4); // this is a ionized h2O with a hole in its last orbital
  
  delete decCh1;
  delete decCh2;
  delete decCh3;
  
  decCh1 = new G4MolecularDecayChannel();
  
  product.SetDefinition( (G4ParticleDefinition*)GetMolecule("OH") );

  product2.SetDefinition( (G4ParticleDefinition*)GetMolecule("H3O") );
  product2.RemoveElectron(5); // I remove the e- of the 5th level, that is the LUMO of H3O+

  decCh1->AddProduct(product);//it should work.
  decCh1->AddProduct(product2);
  decCh1->SetProbability(1);
  
  
  water->AddeConfToExcitedState("ion5", *occ3);
  water->AddDecayChannel("ion5",*decCh1);
  // to this electronic configuration should be associated a decay time of 10e-15 s should the process do it on the dynamic object? the dyn object.
  
  
  //
  // ******************************************
  // *              Building others           *
  // ******************************************
  //
  
  
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAGenericMoleculeManager * G4DNAGenericMoleculeManager::theInstance(0);

