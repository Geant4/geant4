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
/// \file exoticphysics/phonon/src/XPhysicalLattice.cc
/// \brief Implementation of the XPhysicalLattice class
//
// $Id$
//
#include "XPhysicalLattice.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

XPhysicalLattice::XPhysicalLattice(){
  fLattice=NULL;
  fVolume=NULL;
  fTheta=0;
  fPhi=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XPhysicalLattice::XPhysicalLattice(G4VPhysicalVolume* Vol, XLogicalLattice* Lat){
  fLattice=Lat;
  fVolume=Vol;
  fA=fLattice->GetAnhDecConstant();
  fB=fLattice->GetScatteringConstant();
  fDosL=fLattice->GetLDOS();
  fDosST=fLattice->GetSTDOS();
  fDosFT=fLattice->GetFTDOS();
  fBeta=fLattice->GetBeta();
  fGamma=fLattice->GetGamma();
  fLambda=fLattice->GetLambda();
  fMu=fLattice->GetMu();

  G4RotationMatrix *rot = fVolume->GetObjectRotation();

  fGlobalToLocal = G4AffineTransform(*rot);
  fLocalToGlobal = fGlobalToLocal.Invert();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XPhysicalLattice::~XPhysicalLattice(){;}

void XPhysicalLattice::SetDynamicalConstants(double Beta, double Gamma, double Lambda, double Mu)
{
  fBeta=Beta;
  fGamma=Gamma;
  fLambda=Lambda;
  fMu=Mu;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XPhysicalLattice::SetScatteringConstant(G4double a){
  fA=a;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XPhysicalLattice::SetAnhDecConstant(G4double b){
  fB=b;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XPhysicalLattice::SetLDOS(double LDOS){
  fDosL=LDOS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XPhysicalLattice::SetSTDOS(double STDOS)
{
  fDosST = STDOS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XPhysicalLattice::SetFTDOS(double FTDOS){
  fDosFT = FTDOS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XPhysicalLattice::GetBeta(){
  return fBeta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XPhysicalLattice::GetGamma(){
  return fGamma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

double XPhysicalLattice::GetLambda(){
  return fLambda;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XPhysicalLattice::GetMu()
{
  return fMu;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double XPhysicalLattice::GetScatteringConstant()
{
  return fB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double XPhysicalLattice::GetAnhDecConstant()
{
  return fA;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XPhysicalLattice::GetLDOS(){
  return fDosL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XPhysicalLattice::GetSTDOS(){
  return fDosST;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XPhysicalLattice::GetFTDOS()
{
  return fDosFT;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



///////////////////////////////
//Loads the group velocity in m/s
/////////////////////////////
double XPhysicalLattice::MapKtoV(int polarizationState, G4ThreeVector k){
  double groupVelocity;

  k.rotate(G4ThreeVector(0,1,0), fTheta).rotate(G4ThreeVector(0,0,1), fPhi);
  groupVelocity = fLattice->MapKtoV(polarizationState, k);
  k.rotate(G4ThreeVector(0,0,1), -fPhi).rotate(G4ThreeVector(0,1,0), -fTheta);

  return groupVelocity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



///////////////////////////////
//Loads the normalized direction vector along VG
///////////////////////////////
G4ThreeVector XPhysicalLattice::MapKtoVDir(int polarizationState, G4ThreeVector k){

  G4ThreeVector GroupVelocity;

  k=k.rotate(G4ThreeVector(0,1,0), fTheta).rotate(G4ThreeVector(0,0,1), fPhi);
  GroupVelocity = fLattice->MapKtoVDir(polarizationState, k);  

  return GroupVelocity.rotate(G4ThreeVector(0,0,1), -fPhi).rotate(G4ThreeVector(0,1,0), -fTheta).unit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4VPhysicalVolume* XPhysicalLattice::GetVolume(){
  return fVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XPhysicalLattice::SetPhysicalVolume(G4VPhysicalVolume* Vol){
  fVolume=Vol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XPhysicalLattice::SetLatticeOrientation(G4double t_rot, G4double p_rot){
  fTheta=t_rot;
  fPhi= p_rot;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XPhysicalLattice::SetMillerOrientation(int l, int k, int n){
  fTheta=pi/2-std::atan2(n+0.000001,l+0.000001)*rad;
  fPhi=pi/2-std::atan2(l+0.000001,k+0.000001)*rad;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XPhysicalLattice::SetXLogicalLattice(XLogicalLattice* Lat){
  fLattice=Lat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

