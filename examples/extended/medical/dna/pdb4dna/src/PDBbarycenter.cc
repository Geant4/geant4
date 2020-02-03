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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file PDBbarycenter.cc
/// \brief Implementation of the PDBbarycenter class

#include "PDBbarycenter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Barycenter::Barycenter():fBaryNum(0),fRadius(0),
fCenterX(0),fCenterY(0),fCenterZ(0),
fCenterBaseX(0),fCenterBaseY(0),fCenterBaseZ(0),
fCenterSugarX(0),fCenterSugarY(0),fCenterSugarZ(0),
fCenterPhosphateX(0),fCenterPhosphateY(0),fCenterPhosphateZ(0),
fpNext(0)
{
  for ( int i = 0; i < 33; ++i )
  {
    fDistanceTab[i] = 0.; //initialization
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Barycenter::Barycenter(int bNum,double x,double y,double z,
    double Bx,double By, double Bz, //Base barycenter coordinates
    double Sx,double Sy, double Sz, //Sugar barycenter coordinates
    double Px,double Py, double Pz)//Phosphate barycenter coordinates
{
  fBaryNum=bNum;//!<Barycenter number
  fRadius=0;
  fCenterX=x;
  fCenterY=y;
  fCenterZ=z;
  for ( int i = 0; i < 33; ++i )
  {
    fDistanceTab[i] = 0.; //initialization
  }
  fCenterBaseX=Bx;           //!< "X coordinate" of this Base Barycenter
  fCenterBaseY=By;           //!< "Y coordinate" of this Base Barycenter
  fCenterBaseZ=Bz;           //!< "Z coordinate" of this Base Barycenter
  fCenterSugarX=Sx;          //!< "X coordinate" of this Sugar Barycenter
  fCenterSugarY=Sy;          //!< "Y coordinate" of this Sugar Barycenter
  fCenterSugarZ=Sz;          //!< "Z coordinate" of this Sugar Barycenter
  fCenterPhosphateX=Px;      //!< "X coordinate" of this Phosphate Barycenter
  fCenterPhosphateY=Py;      //!< "Y coordinate" of this Phosphate Barycenter
  fCenterPhosphateZ=Pz;      //!< "Z coordinate" of this Phosphate Barycenter
  fpNext=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Barycenter *Barycenter::GetNext()
{
  return fpNext;
}

int Barycenter::GetID()
{
  return fBaryNum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Barycenter::SetNext(Barycenter *barycenterNext)
{
  fpNext=barycenterNext;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Barycenter::SetDistance(int i, double dist)
{
  fDistanceTab[i]=dist;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Barycenter::GetDistance(int i)
{
  return fDistanceTab[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Barycenter::SetRadius(double rds)
{
  fRadius=rds;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Barycenter::GetRadius()
{
  return fRadius;
}

