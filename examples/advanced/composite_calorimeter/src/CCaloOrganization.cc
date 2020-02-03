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
///////////////////////////////////////////////////////////////////////////////
// File: CCaloOrganization.cc
// Description: Packing, unpacking and other related utilities for 
//              calorimetric numbering schema
///////////////////////////////////////////////////////////////////////////////
#include "CCaloOrganization.hh"

//#define debug

unsigned int CCaloOrganization::packindex(G4int det, G4int z, G4int eta, 
                                          G4int phi) const {
  //So this is the actual encoding of the index:
  //top 4 bits encode Detector type
  //
  // Should work for all calorimeter with no depth information

  unsigned int idx=(det&15)<<28;   //bits 28-31   (21-23 are free for now)
  idx+=(((z+1)/2)&1)<<20;          //bits 20
  idx+=(eta&1023)<<10;             //bits 10-19
  idx+=(phi&1023);                 //bits  0-9
#ifdef debug
  G4cout << " ECAL packing " << det << " " << z << " " << eta << " " << phi 
         << "  into " << idx << G4endl;
#endif
  return idx;
}

unsigned int CCaloOrganization::packindex(G4int det, G4int depth, G4int z, G4int eta, 
                                          G4int phi) const {
  //So this is the actual encoding of the index:
  //top 4 bits encode Detector type
  //next 4 bits encode depth information
  // Should work for all calorimeter with no depth information

  unsigned int idx=(det&15)<<28;  //bits 28-31   (21-23 are free for now)
  idx+=(depth&15)<<24;            //bits 24-27  
  idx+=(z&1)<<20;                 //bits 20
  idx+=(eta&1023)<<10;            //bits 10-19
  idx+=(phi&1023);                //bits  0-9
#ifdef debug
  G4cout << " HCAL packing " << det << " " << depth << " " << z << " " << eta 
         << " " << phi  << "  into " << idx << G4endl;
#endif
  return idx;
}


void CCaloOrganization::unpackindex(const unsigned int& idx, G4int& det, G4int& z, 
                                    G4int& eta, G4int& phi) const {
  det = (idx>>28)&15;
  z   = (idx>>20)&1;
  z   = 2*z-1;
  eta = (idx>>10)&1023;
  phi = (idx&1023);
}


void CCaloOrganization::unpackindex(const unsigned int& idx, G4int& det, 
                                    G4int& depth, G4int& z, G4int& eta, 
                                    G4int& phi) const {
  det = (idx>>28)&15;
  depth=(idx>>24)&15;
  z   = (idx>>20)&1;
  eta = (idx>>10)&1023;
  phi = (idx&1023);
}


G4int CCaloOrganization::getUnitWithMaxEnergy(std::map<G4int,G4float,std::less<G4int> >& themap){

  //look for max
  G4int UnitWithMaxEnergy = 0;
  G4float maxEnergy = 0.;
        
  for(std::map<G4int,G4float,std::less<G4int> >::iterator iter = themap.begin();
      iter != themap.end(); iter++){
            
    if(        maxEnergy < (*iter).second) {
      maxEnergy = (*iter).second;        
      UnitWithMaxEnergy = (*iter).first;
    }                                
  }        
#ifdef debug
  G4cout << " *** max energy of " << maxEnergy << " MeV was found in Unit id "
         << UnitWithMaxEnergy;
  G4int det,z,eta,phi;
  unpackindex(UnitWithMaxEnergy, det, z, eta, phi);
  G4cout << " corresponding to z= " << z << " eta= " << eta << " phi = " << phi
         << G4endl;
#endif
  return UnitWithMaxEnergy;

}


G4float CCaloOrganization::energyInMatrix(G4int nCellInEta, G4int nCellInPhi,
                                        G4int crystalWithMaxEnergy, 
                                        std::map<G4int,G4float,std::less<G4int> >& themap){

  G4int det,z,eta,phi;
  this->unpackindex(crystalWithMaxEnergy, det, z, eta, phi);
  G4int ncristals=0;
        
  G4int goBackInEta = nCellInEta/2;
  G4int goBackInPhi = nCellInPhi/2;
  G4int startEta = eta-goBackInEta;
  G4int startPhi = phi-goBackInPhi;

  G4float totalEnergy = 0.;
  
  for(G4int ieta=startEta; ieta<startEta+nCellInEta; ieta++){
    for(G4int iphi=startPhi; iphi<startPhi+nCellInPhi; iphi++){
      
      G4int index = this->packindex(det,z,ieta,iphi);
      totalEnergy += themap[index];
      ncristals+=1;
#ifdef debug
      G4cout << "ieta - iphi - E = " << ieta << "  " << iphi << " " 
             << themap[index] << G4endl;
#endif
    }
  }
        
#ifdef debug
  G4cout << "Energy in " << nCellInEta << " cells in eta times "
         << nCellInPhi << " cells in phi matrix = " << totalEnergy
         << " for " << ncristals << " crystals" << G4endl;
#endif
  return totalEnergy;

}   
