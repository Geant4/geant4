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
// $Id: $
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class implementation file
//
// File name:     G4PWATotalXsecTable
//
// Author:        Mihaly Novak
//
// Creation date: 18.05.2015
//
// Class description:
//   Class to load and handle elastic, first and second transport cross sections
//   precomputed by using ELSEPA [1] in the 100 eV - 1 GeV kinetic and Z = 1-103
//   energy range for electrons and positrons.G4PWATotalXsecZ is responsible to 
//   to handle cross sections by individual Z that are used in the current 
//   geometry and G4PWATotalXsecTable is a collection of G4PWATotalXsecZ objects. 
//
// Modifications:
//
// References:
//   [1] Francesc Salvat, Aleksander Jablonski, Cedric J Powell,
//       ELSEPA—Dirac partial-wave calculation of elastic scattering of electrons 
//       and positrons by atoms, positive ions and molecules,
//       Computer physics communications; 165, 2, (2005)
//
// -----------------------------------------------------------------------------


#include "G4PWATotalXsecTable.hh"


#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

////////////////////////////////////////////////////////////////////////////////
//  G4PWATotalXsecZ: sub-class for PWA xsec data that belong to a given Z number
//////////////////////////////////////////////////////////////////////////////// 

////////////////////////////////////////////////////////////////////////////////
const G4double G4PWATotalXsecZ::fgPWATotalXsecEnergyGrid[]={
    // energy bin values for total elastic, first and second transport cross scetions in MeV
    1.00000000e-04, 1.16591440e-04, 1.35935639e-04, 1.58489319e-04, 1.84784980e-04, 2.15443469e-04, 2.51188643e-04, 2.92864456e-04, 
    3.41454887e-04, 3.98107171e-04, 4.64158883e-04, 5.41169527e-04, 6.30957344e-04, 7.35642254e-04, 8.57695899e-04, 1.00000000e-03, 
    1.16591440e-03, 1.35935639e-03, 1.58489319e-03, 1.84784980e-03, 2.15443469e-03, 2.51188643e-03, 2.92864456e-03, 3.41454887e-03, 
    3.98107171e-03, 4.64158883e-03, 5.41169527e-03, 6.30957344e-03, 7.35642254e-03, 8.57695899e-03, 1.00000000e-02, 1.16591440e-02, 
    1.35935639e-02, 1.58489319e-02, 1.84784980e-02, 2.15443469e-02, 2.51188643e-02, 2.92864456e-02, 3.41454887e-02, 3.98107171e-02, 
    4.64158883e-02, 5.41169527e-02, 6.30957344e-02, 7.35642254e-02, 8.57695899e-02, 1.00000000e-01, 1.16591440e-01, 1.35935639e-01, 
    1.58489319e-01, 1.84784980e-01, 2.15443469e-01, 2.51188643e-01, 2.92864456e-01, 3.41454887e-01, 3.98107171e-01, 4.64158883e-01, 
    5.41169527e-01, 6.30957344e-01, 7.35642254e-01, 8.57695899e-01, 1.00000000e+00, 1.16591440e+00, 1.35935639e+00, 1.58489319e+00, 
    1.84784980e+00, 2.15443469e+00, 2.51188643e+00, 2.92864456e+00, 3.41454887e+00, 3.98107171e+00, 4.64158883e+00, 5.41169527e+00, 
    6.30957344e+00, 7.35642254e+00, 8.57695899e+00, 1.00000000e+01, 1.16591440e+01, 1.35935639e+01, 1.58489319e+01, 1.84784980e+01, 
    2.15443469e+01, 2.51188643e+01, 2.92864456e+01, 3.41454887e+01, 3.98107171e+01, 4.64158883e+01, 5.41169527e+01, 6.30957344e+01, 
    7.35642254e+01, 8.57695899e+01, 1.00000000e+02, 1.16591440e+02, 1.35935639e+02, 1.58489319e+02, 1.84784980e+02, 2.15443469e+02, 
    2.51188643e+02, 2.92864456e+02, 3.41454887e+02, 3.98107171e+02, 4.64158883e+02, 5.41169527e+02, 6.30957344e+02, 7.35642254e+02, 
    8.57695899e+02, 1.00000000e+03
 };

////////////////////////////////////////////////////////////////////////////////
G4PWATotalXsecZ::G4PWATotalXsecZ(G4int Z){
  G4int nn = fgNumTotalXsecBins*6;
  for(G4int i=0; i<nn; ++i) {
    fPWAXsecs[i] = 0.0;
    fInterpParamA[i] = 0.0;
    fInterpParamB[i] = 0.0;
  }
  LoadPWATotalXsecZ(Z); 
}

////////////////////////////////////////////////////////////////////////////////
void G4PWATotalXsecZ::LoadPWATotalXsecZ(G4int Z){
 G4double dum;
 char fname[512];
 char* path = getenv("G4LEDATA");
 if (!path) {
    G4Exception("G4PWATotalXsecZ::LoadPWATotalXsecZ()","em0006",
		FatalException,
		"Environment variable G4LEDATA not defined");
    return;
 }

 std::string pathString(path);
 sprintf(fname,"%s/msc_GS/xsecs/xsecs_%d",path,Z);
 std::ifstream infile(fname,std::ios::in);
 if(!infile.is_open()){
    char msgc[512];
    sprintf(msgc,"  Total PWA xsection %s not found.",fname);
    G4Exception("G4PWATotalXsecZ::LoadPWATotalXsecZ()","em0006",
		FatalException,
		msgc);
    return;
 }
 G4double dummy;

 for(G4int i=0; i<fgNumTotalXsecBins; ++i)
   for(G4int j=0; j<7; ++j)
     if(j==0) infile >> dum;
     else {
       // load pwa xsection that are stored in cm2 units in file and change to 
       // Geant4 internal length2 units 
       infile >> dummy;
       fPWAXsecs[(j-1)*fgNumTotalXsecBins+i] = dummy*CLHEP::cm2;        
     }
 infile.close();

 // compute log-log linear intrp. parameters
 for(G4int i=0; i<fgNumTotalXsecBins-1; ++i)
    for(G4int k=0; k<6; ++k) {
       G4int j = k*fgNumTotalXsecBins+i; 
       G4double val2 = fPWAXsecs[j+1];
       G4double val1 = fPWAXsecs[j];

       fInterpParamA[j] =  G4Log(val2/val1)/G4Log(fgPWATotalXsecEnergyGrid[i+1]/fgPWATotalXsecEnergyGrid[i]);
       fInterpParamB[j] =  G4Exp(G4Log(val1) - fInterpParamA[j]*G4Log(fgPWATotalXsecEnergyGrid[i]));
 }
}

////////////////////////////////////////////////////////////////////////////////
// Get the index of the lower energy bin edge
G4int G4PWATotalXsecZ::GetPWATotalXsecEnergyBinIndex(G4double energy) const {
  // log(fgPWATotalXsecEnergyGrid[0]);
  const G4double lne0 = -9.21034037197618e+00; 
  // 1./log(fgPWATotalXsecEnergyGrid[i+1]/fgPWATotalXsecEnergyGrid[i]); 
  const G4double invlnde =  6.51441722854880e+00;   

  return (G4int)((G4Log(energy)-lne0)*invlnde);
}

////////////////////////////////////////////////////////////////////////////////
//  j-dependent type interploated cross section in Geant4 internal length2 unit
// energy is assumed to be in [MeV]
G4double G4PWATotalXsecZ::GetInterpXsec(G4double energy, G4int elowindex, G4int j) const { 
   // protection : out of energy grid range 
   if(energy < GetLowestEnergy())
      return GetLowestXsecValue(j);
   if(energy >= GetHighestEnergy())
      return GetHighestXsecValue(j);     

   // normal case log-log linear intrp.
  G4int k = j*fgNumTotalXsecBins+elowindex; 
  return G4Exp(G4Log(energy)*fInterpParamA[k])*fInterpParamB[k];
}

////////////////////////////////////////////////////////////////////////////////
G4double G4PWATotalXsecZ::GetInterpXsec(G4double energy, G4int j) const {
   // protection : out of energy grid range 
   if(energy < GetLowestEnergy())
      return GetLowestXsecValue(j);
   if(energy >= GetHighestEnergy())
      return GetHighestXsecValue(j);     

   // normal case log-log linear intrp.
  G4int elowindex = GetPWATotalXsecEnergyBinIndex(energy); 
  G4int k = j*fgNumTotalXsecBins+elowindex; 
  return G4Exp(G4Log(energy)*fInterpParamA[k])*fInterpParamB[k];
}


////////////////////////////////////////////////////////////////////////////////
//  G4PWATotalXsecTable
////////////////////////////////////////////////////////////////////////////////
 
G4PWATotalXsecZ* G4PWATotalXsecTable::fgPWATotalXsecTable[fgNumZet] = {0};

////////////////////////////////////////////////////////////////////////////////
G4PWATotalXsecTable::~G4PWATotalXsecTable(){
     for(G4int i = 0; i < fgNumZet; ++i) 
       if(fgPWATotalXsecTable[i]) {
         delete fgPWATotalXsecTable[i];
         fgPWATotalXsecTable[i] = 0; 
       } 
}

////////////////////////////////////////////////////////////////////////////////
void G4PWATotalXsecTable::Initialise(){
   G4int isUsedZ[fgNumZet] ={0};  //xsec data available up to fgNumZet Z-number
   // check used elements
   G4MaterialTable *theMaterialTable = G4Material::GetMaterialTable();
   for(unsigned int imat = 0; imat < theMaterialTable->size(); ++imat) {
     const G4ElementVector  *theElemVect = ((*theMaterialTable)[imat])->GetElementVector(); 
     for(unsigned int ielem = 0; ielem < theElemVect->size(); ++ielem) {
       G4int zet = G4lrint((*theElemVect)[ielem]->GetZ());
       zet = zet>fgNumZet ? fgNumZet : zet;  
       if(!isUsedZ[zet-1])
         isUsedZ[zet-1] = 1;
     }
   }

   for(G4int i = 0; i < fgNumZet; ++i)
     if(isUsedZ[i] && !fgPWATotalXsecTable[i])       // used but not there yet -> load it
       fgPWATotalXsecTable[i] = new G4PWATotalXsecZ(i+1); 
     else if(!isUsedZ[i] && fgPWATotalXsecTable[i]) { // there but not used now -> delete
       delete fgPWATotalXsecTable[i];
       fgPWATotalXsecTable[i] = 0;  
     }
}

