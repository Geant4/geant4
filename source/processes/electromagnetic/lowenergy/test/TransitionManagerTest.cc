// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      Geant4 Collaboration
//
//      File name:     TransitionManagerTest.cc
//
//      Authors:        Alfonso Mantero (alfonso.mantero@ge.infn.it)
//                      Elena Guardincerri (elena.guardincerri@ge.infn.it) 
//
//      Creation date: 1/05/2001
// -------------------------------------------------------------------


#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4AtomicShell.hh"
#include "G4AtomicTransitionManager.hh"

int main()
{
   //--------- Materials definition ---------

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
  G4int nElements = theMaterialTable->length();

  static G4AtomicTransitionManager* manager = G4AtomicTransitionManager::Instance();


 G4int z = 0;
 
 G4int shellIdentifier = 0;
 
 G4cout<<"Insert Z"<<G4endl;
 
 G4cin>> z;

 
 G4cout<<"Insert the starting shell"<<G4endl;
 
 G4cin >> shellIdentifier;

 const G4AtomicShell* theAtomicShell = manager->Shell(z,shellIdentifier);


 G4int theShellId = theAtomicShell->ShellId();

 G4cout <<"theShellId ="<< theShellId <<G4endl;

G4double theBindingEnergy = theAtomicShell->BindingEnergy();

 G4cout<<"theBindingEnergy = "<< theBindingEnergy << G4endl;

 const G4DataVector theTransitionProbabilities = theAtomicShell->TransitionProbabilities();

 G4int sizeProb=theTransitionProbabilities.size();

 for (G4int i = 0;i < sizeProb; i++){

   G4cout << "theTransitionProbabilities ["<<i<<"] = " << theTransitionProbabilities[i] << G4endl;

 }

 const G4DataVector theTransitionEnergies = theAtomicShell->TransitionEnergies();

G4int sizeEn = theTransitionEnergies.size();

 for (G4int j = 0;j < sizeEn; j++){

 G4cout << "theTransitionEnergies ["<<j<<"] = " << theTransitionEnergies[j] << G4endl;

 }

const G4DataVector theTransFinalShellIdentifiers = theAtomicShell->TransFinalShellIdentifiers();

G4int sizeId = theTransFinalShellIdentifiers.size();

 for (G4int k = 0;k < sizeId; k++){ 

 G4cout << "theTransFinalShellIdentifiers ["<<k<<"] = " <<theTransFinalShellIdentifiers[k] << G4endl;
 
 }

 G4int numberOfShells = manager->NumberOfShells(z);

 G4cout<<" The number of shells is :"<<numberOfShells<<G4endl;

 G4double totalRadiativeTransitionProbability = manager->TotalRadiativeTransitionProbability(z,shellIdentifier);

 G4cout<<" The total radiative transition probability is: "<<totalRadiativeTransitionProbability<< G4endl;

 G4double totalNonRadiativeTransitionProbability = manager->TotalNonRadiativeTransitionProbability(z,shellIdentifier);

G4cout<<" The total non radiative transition probability is: "<<totalNonRadiativeTransitionProbability<< G4endl;

 G4cout << "END OF THE MAIN PROGRAM" << G4endl;
 
}
