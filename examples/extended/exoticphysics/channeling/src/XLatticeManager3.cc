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

#include "XLatticeManager3.hh"
#include "G4VPhysicalVolume.hh"

//int XLatticeManager3::fTotalLattices = 0;
XLatticeManager3* XLatticeManager3::LM;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XLatticeManager3::XLatticeManager3()
{
  fTotalLattices = 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLatticeManager3::~XLatticeManager3()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XLatticeManager3* XLatticeManager3::GetXLatticeManager(){

  //if no lattice manager exists, create one.
  if(!LM) LM = new XLatticeManager3();

  //return pointer to single existing lattice manager
  return LM;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

bool XLatticeManager3::RegisterLattice(XPhysicalLattice* Lat){
  if(fTotalLattices<MAXLAT){
    fTotalLattices++;
    //using "fTotalLattices-1" so that first lattice corresponds to index 0
    fLatticeList[fTotalLattices-1]=Lat;  
    G4cout<<"\nXLatticeManager3::registerLattice: Registering Lattice.";
    G4cout<< "Total number of lattices:"<<fTotalLattices<<"\n"<<endl;

    return true; 
  }
  
  G4cout<<"\nXLatticeManager::RegisterLattice(XPhysicalLattice*):";
  G4cout<<"Maximum number of lattices MAXLAT exceeded."<<endl;
 
  return false;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicalLattice* XLatticeManager3::GetXPhysicalLattice(G4VPhysicalVolume* Vol){
//returns a pointer to the PhysicalLattice associated with Vol

    for(int counter=0;counter<fTotalLattices;counter++){
      if(fLatticeList[counter]->GetVolume()==Vol) {
        return fLatticeList[counter]; //found matching lattice
      }      
    }
    return fLatticeList[0];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

bool XLatticeManager3::HasLattice(G4VPhysicalVolume* Vol){
  //return true if Vol has a physical lattice

    for(int counter=0;counter<fTotalLattices;counter++){
      if(fLatticeList[counter]->GetVolume()==Vol) {
        return true; //found matching lattice
      }      
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

double XLatticeManager3::MapKtoV(G4VPhysicalVolume* Vol,
                                 int polarizationState,
                                 const G4ThreeVector & k)
{
  //Given the phonon wave vector k, phonon physical volume Vol 
  //and polarizationState(0=LON, 1=FT, 2=ST), 
  //returns phonon velocity in m/s

  if((Vol==NULL)&&(fTotalLattices>0))
      return fLatticeList[0]->MapKtoV(polarizationState, k);
  for(int counter=0;counter<fTotalLattices;counter++){
    if(fLatticeList[counter]->GetVolume()==Vol) {
      return fLatticeList[counter]->MapKtoV(polarizationState, k);
    }
  }
  G4cout<<"\nXLatticeManager::MapKtoV: Found no matching lattices for "
    <<Vol->GetName()<<". Total number of lattices is "<<fTotalLattices<<endl;
return 300;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XLatticeManager3::MapKtoVDir(G4VPhysicalVolume* Vol,
                                           int polarizationState,
                                           const G4ThreeVector & k)
{
  //Given the phonon wave vector k, phonon physical volume Vol 
  //and polarizationState(0=LON, 1=FT, 2=ST), 
  //returns phonon propagation direction as dimensionless unit vector

  if((Vol==NULL)&&(fTotalLattices>0))
      return fLatticeList[0]->MapKtoVDir(polarizationState, k);
  for(int counter=0;counter<fTotalLattices;counter++){
    if(fLatticeList[counter]->GetVolume()==Vol) {
      if(counter!=0)
          G4cout<<"\nLattiveManager2::MapKtoV:";
          G4cout<<"returning group velocity from lattise position: "<<counter;
      return fLatticeList[counter]->MapKtoVDir(polarizationState, k);
    }
  }
return G4ThreeVector(1,0,0);
}

