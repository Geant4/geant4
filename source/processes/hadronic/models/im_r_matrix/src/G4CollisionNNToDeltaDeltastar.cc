//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4CollisionNNToDeltaDeltastar.cc,v 1.1 2003/10/07 12:37:35 hpw Exp $ //

#include "globals.hh"
#include "G4CollisionNNToDeltaDeltastar.hh"
#include "G4KineticTrack.hh"
#include "G4VCrossSectionSource.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4XAqmElastic.hh"
#include "G4AngularDistribution.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4KineticTrackVector.hh"
#include "G4ParticleTable.hh"
#include "G4CollisionVector.hh"
#include "G4CollisionNNToDeltaDeltastar.hh"
#include "G4CollisionNNToDeltaDelta1600.hh"
#include "G4CollisionNNToDeltaDelta1620.hh"
#include "G4CollisionNNToDeltaDelta1700.hh"
#include "G4CollisionNNToDeltaDelta1900.hh"
#include "G4CollisionNNToDeltaDelta1905.hh"
#include "G4CollisionNNToDeltaDelta1910.hh"
#include "G4CollisionNNToDeltaDelta1920.hh"
#include "G4CollisionNNToDeltaDelta1930.hh"
#include "G4CollisionNNToDeltaDelta1950.hh"

G4CollisionNNToDeltaDeltastar::G4CollisionNNToDeltaDeltastar()
{ 
  
  // 1600
  G4CollisionComposite::AddComponent(new G4CollisionNNToDeltaDelta1600()); 
  //1620 
  G4CollisionComposite::AddComponent(new G4CollisionNNToDeltaDelta1620()); 
  //1700 
  G4CollisionComposite::AddComponent(new G4CollisionNNToDeltaDelta1700()); 
  //1900 
  G4CollisionComposite::AddComponent(new G4CollisionNNToDeltaDelta1900()); 
  //1905 
  G4CollisionComposite::AddComponent(new G4CollisionNNToDeltaDelta1905()); 
  //1910 
  G4CollisionComposite::AddComponent(new G4CollisionNNToDeltaDelta1910()); 
  //1920 
  G4CollisionComposite::AddComponent(new G4CollisionNNToDeltaDelta1920()); 
  //1930 
  G4CollisionComposite::AddComponent(new G4CollisionNNToDeltaDelta1930()); 
  //1950 
  G4CollisionComposite::AddComponent(new G4CollisionNNToDeltaDelta1950()); 

}

