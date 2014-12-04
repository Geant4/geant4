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
#include "G4LEPTSDistribution.hh"


G4LEPTSDistribution::G4LEPTSDistribution()
{
}


void G4LEPTSDistribution::ReadFile(G4String fileName) {

  G4int eB, out, out2;
  float  float_data1,float_data2;
  G4double sum, esum;
  FILE * fp;

  for (eB=0;eB<10000;eB++){
    E[eB]=0.0;
    f[eB]=0.0;
    F[eB]=0.0;
    eF[eB]=0.0;
  }

  if ((fp=fopen(fileName.c_str(), "r"))==NULL){
    //G4cout << "Error reading " << fileName << G4endl;
    NoBins = 0;
    bFileFound = false;
    return;
  }
  else{
    bFileFound = true;
    //    G4cout << "Read Distro (" << fileName << ") " << G4endl;
    out=1;
    eB=1;
    while (out==1){
      out  = fscanf(fp,"%f \n",&float_data1);
      out2 = fscanf(fp,"%f \n",&float_data2);
      if (out==1 && out2==1){
	E[eB]=(G4double)float_data1;
	f[eB]=(G4double)float_data2;
	eB++;
      }
    }

    fclose(fp);
  }

  NoBins=eB-1;  //=1272+1 or 9607+1;

  if( NoBins >= NMAX )
    printf("ERROR !!!!  Eloss NoBins= %d \n", NoBins);

  sum=0.0;
  esum=0.0;
  for (eB=0;eB<=NoBins;eB++) {
    if( f[eB] > 0) {
      sum+=f[eB];
      esum+=E[eB]*f[eB];
    }
    F[eB]=sum;
    eF[eB]=esum;
  }

  //  if( verboseLevel >= 1 ) G4cout << "Norm: " << F[NoBins] << " NoBins: "<< NoBins << G4endl;
  
  for (eB=0;eB<=NoBins;eB++) {
    eF[eB] = eF[eB]/F[eB];
    F[eB] = F[eB]/F[NoBins];
  }
  //for (eB=0;eB<=NoBins;eB++)
  //G4cout << "eff " << E[eB] << " " << f[eB] << " " << F[eB] << "\n";
}

G4bool G4LEPTSDistribution::ReadFile( FILE* fp, G4int nData ) 
{

  G4int eB, out, out2;
  float  float_data1,float_data2;
  G4double sum, esum;

  for (eB=0;eB<10000;eB++){
    E[eB]=0.0;
    f[eB]=0.0;
    F[eB]=0.0;
    eF[eB]=0.0;
  }

  bFileFound = true;
  out=1;
  eB=1;

  for( G4int id = 0; id < nData; id++ ){	  
    out  = fscanf(fp,"%f \n",&float_data1);
    out2 = fscanf(fp,"%f \n",&float_data2);
    if (out==1 && out2==1){
      E[eB]=(G4double)float_data1;
      f[eB]=(G4double)float_data2;
      eB++;
    }else{
      return 1;
    }
  }
    
  NoBins=eB-1;  //=1272+1 or 9607+1;

  if( NoBins >= NMAX )
    printf("ERROR !!!!  Eloss NoBins= %d \n", NoBins);

  sum=0.0;
  esum=0.0;
  for (eB=0;eB<=NoBins;eB++) {
    if( f[eB] > 0) {
      sum+=f[eB];
      esum+=E[eB]*f[eB];
    }
    F[eB]=sum;
    eF[eB]=esum;
  }

  //if( verboseLevel >= 1 ) G4cout << "Norm: " << F[NoBins] << " NoBins: "<< NoBins << G4endl;
  
  for (eB=0;eB<=NoBins;eB++) {
    eF[eB] = eF[eB]/F[eB];
    F[eB] = F[eB]/F[NoBins];
  }
  //for (eB=0;eB<=NoBins;eB++)
  //G4cout << "eff " << E[eB] << " " << f[eB] << " " << F[eB] << "\n";

  return 0;
}



G4double G4LEPTSDistribution::Sample( G4double eMin, G4double eMax) {
// Sample Energy from Cumulative distr. G4interval [eMin, eMax]

  if( eMin > eMax) return 0.0;

  G4int i,j,k=0, iMin, iMax;

  i=0; j=NoBins;
  while ((j-i)>1) {
    k=(i+j)/2;
    if( E[k] < eMax ) i=k;
    else              j=k;
  }
  iMax = i;

  i=0; j=NoBins;
  while ((j-i)>1) {
    k=(i+j)/2;
    if( E[k] < eMin ) i=k;
    else              j=k;
  }
  iMin = i;

  G4double rnd = F[iMin] + (F[iMax] - F[iMin]) * G4UniformRand();

  i=0; j=NoBins;
  while ((j-i)>1) {
    k=(i+j)/2;
    if( F[k]<rnd) i=k;
    else          j=k;
  }

  G4double Sampled = E[k];

  if(      Sampled < eMin) Sampled = eMin;
  else if( Sampled > eMax) Sampled = eMax;

  return Sampled;
}
