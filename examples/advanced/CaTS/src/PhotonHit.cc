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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file PhotonHit.cc
/// \brief Implementation of the CaTS::PhotonHit class

// Geant4 headers
#include "PhotonHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <G4Point3D.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>
// project headers
#include "PhotonHit.hh"

template <class Type>
class G4Allocator;
G4ThreadLocal G4Allocator<PhotonHit>* PhotonHitAllocator = nullptr;

PhotonHit::PhotonHit()
  : G4VHit()
{}

PhotonHit::PhotonHit(unsigned iid, unsigned ipid, G4double iwavelength,
                     G4double itime, G4ThreeVector iposition,
                     G4ThreeVector idirection, G4ThreeVector ipolarization)
  : G4VHit()
{
  fid           = iid;
  fpid          = ipid;
  fwavelength   = iwavelength;
  ftime         = itime;
  fposition     = iposition;
  fdirection    = idirection;
  fpolarization = ipolarization;
}

PhotonHit::PhotonHit(const PhotonHit& right)
  : G4VHit()
{
  fid           = right.fid;
  fpid          = right.fpid;
  fwavelength   = right.fwavelength;
  ftime         = right.ftime;
  fposition     = right.fposition;
  fdirection    = right.fdirection;
  fpolarization = right.fpolarization;
}

const PhotonHit& PhotonHit::operator=(const PhotonHit& right)
{
  fid           = right.fid;
  fpid          = right.fpid;
  fwavelength   = right.fwavelength;
  ftime         = right.ftime;
  fposition     = right.fposition;
  fdirection    = right.fdirection;
  fpolarization = right.fpolarization;
  return *this;
}

G4bool PhotonHit::operator==(const PhotonHit& right) const
{
  return (this == &right) ? true : false;
}

void PhotonHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fposition);
    circle.SetScreenSize(2.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1., 0., 0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}
