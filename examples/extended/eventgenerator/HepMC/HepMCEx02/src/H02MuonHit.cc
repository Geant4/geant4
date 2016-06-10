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
/// \file eventgenerator/HepMC/HepMCEx02/src/H02MuonHit.cc
/// \brief Implementation of the H02MuonHit class
//
//   $Id: H02MuonHit.cc 77801 2013-11-28 13:33:20Z gcosmo $
//
#include <iomanip>
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "H02MuonHit.hh"

G4Allocator<H02MuonHit> H02MuonHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02MuonHit::H02MuonHit()
  : moduleID(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02MuonHit::H02MuonHit(G4int imod, G4String aname,
                     const G4ThreeVector& pxyz,
                     const G4ThreeVector& xyz, G4double atof)
  : moduleID(imod), pname(aname), momentum(pxyz),
    position(xyz), tof(atof)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02MuonHit::~H02MuonHit()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02MuonHit::H02MuonHit(const H02MuonHit& right)
  : G4VHit()
{
  *this= right;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const H02MuonHit& H02MuonHit::operator=(const H02MuonHit& right)
{
  moduleID= right.moduleID;
  pname= right.pname;
  momentum= right.momentum;
  position= right.position;
  tof= right.tof;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int H02MuonHit::operator==(const H02MuonHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02MuonHit::Draw()
{
  const G4double pt_min=20.*GeV;

  G4VVisManager* pVVisManager= G4VVisManager::GetConcreteInstance();
  if(pVVisManager) {
    G4Circle circle(position);
    circle.SetScreenSize(5.);
    circle.SetFillStyle(G4Circle::filled);

    G4Color color, goodColor(1.,0.,0.), badColor(0.,0.,1.);
    if(momentum.perp()>pt_min) color=goodColor;
    else color=badColor;

    G4VisAttributes attribs(color);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02MuonHit::Print()
{
  G4int id= moduleID;
  G4String tag="B";
  if(moduleID >=10) {
    id -=10;
    tag="E";
  }
  G4cout << tag << id << " :" << std::setw(12) << pname.c_str()
         << " : pT=" << std::setprecision(3)  << momentum.perp()/GeV
         << " : TOF=" << std::setprecision(3) << tof/ns
         << " : x="  << std::setprecision(3) << position*(1./m)
         << G4endl;
}
