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
// -------------------------------------------------------------------
//
//      Author:        E.Mendoza
// 
//      Creation date: May 2024
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//  NuDEX code (https://doi.org/10.1016/j.nima.2022.167894)
// 



#include "G4NuDEXRandom.hh"

#if COMPILATIONTYPE == 1
//==============================================================================
G4NuDEXRandom::G4NuDEXRandom(unsigned int seed){
  theRandom=new TRandom2(seed);
}
G4NuDEXRandom::~G4NuDEXRandom(){
  delete theRandom;
}
void G4NuDEXRandom::SetSeed(unsigned int seed){
  theRandom->SetSeed(seed);
}
unsigned int G4NuDEXRandom::GetSeed(){
  return theRandom->GetSeed();
}
G4double G4NuDEXRandom::Uniform(G4double Xmin,G4double Xmax){
  return theRandom->Uniform(Xmin,Xmax);
}
unsigned int G4NuDEXRandom::Integer(unsigned int IntegerMax){
  return theRandom->Integer(IntegerMax);
}
G4double G4NuDEXRandom::Exp(G4double tau){
  return theRandom->Exp(tau);
}
G4double G4NuDEXRandom::Gaus(G4double mean,G4double sigma){
  return theRandom->Gaus(mean,sigma);
}
G4int G4NuDEXRandom::Poisson(G4double mean){
  return theRandom->Poisson(mean);
}
//==============================================================================
void NuDEXException(const char* originOfException, const char* exceptionCode,const char* ){
  std::cout<<" ############## Error in "<<originOfException<<", line "<<exceptionCode<<" ##############"<<std::endl; exit(1);
}
//==============================================================================

#elif COMPILATIONTYPE == 2
//==============================================================================
G4NuDEXRandom::G4NuDEXRandom(unsigned int seed){
  theEngine=new CLHEP::HepJamesRandom(seed);
  theRandFlat=new CLHEP::RandFlat(theEngine);
  theRandExponential=new CLHEP::RandExponential(theEngine);
  theRandGauss=new CLHEP::RandGauss(theEngine);
  theRandPoisson=new CLHEP::RandPoisson(theEngine);
}
G4NuDEXRandom::~G4NuDEXRandom(){

  //delete theRandFlat;
  //delete theRandExponential;
  //delete theRandGauss;
  //delete theRandPoisson;
  //delete theEngine;

}
void G4NuDEXRandom::SetSeed(unsigned int seed){
  theEngine->setSeed(seed);
  theRandGauss->setF(false);
}
unsigned int G4NuDEXRandom::GetSeed(){
  return (unsigned int)theEngine->getSeed();
}
G4double G4NuDEXRandom::Uniform(G4double Xmin,G4double Xmax){
  return theRandFlat->fire(Xmin,Xmax);
}
unsigned int G4NuDEXRandom::Integer(unsigned int IntegerMax){
  return (unsigned int)theRandFlat->fireInt(IntegerMax); //bikerful!!!
}
G4double G4NuDEXRandom::Exp(G4double tau){
  return theRandExponential->fire(tau);
}
G4double G4NuDEXRandom::Gaus(G4double mean,G4double sigma){
  return theRandGauss->fire(mean,sigma);
}
G4long G4NuDEXRandom::Poisson(G4double mean){
  return theRandPoisson->fire(mean);
}
//==============================================================================
void NuDEXException(const char* originOfException, const char* exceptionCode,const char* description){
  G4Exception(originOfException,exceptionCode,FatalException,description);
}
//==============================================================================
#endif
