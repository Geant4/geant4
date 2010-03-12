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
// $Id: G4CollisionNStarNToNN.cc,v 1.3 2010-03-12 15:45:18 gunter Exp $ //

#include "globals.hh"
#include "G4CollisionNStarNToNN.hh"
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
#include "G4CollisionNStarNToNN.hh"
#include "G4ConcreteNStarNToNN.hh"

typedef G4ConcreteNStarNToNN channelType;

G4CollisionNStarNToNN::G4CollisionNStarNToNN()
{ 
  MakeNNStarToNN<channelType, N1400pPC, N1400nPC>::Make(this);
  MakeNNStarToNN<channelType, N1520pPC, N1520nPC>::Make(this);
  MakeNNStarToNN<channelType, N1535pPC, N1535nPC>::Make(this);
  MakeNNStarToNN<channelType, N1650pPC, N1650nPC>::Make(this);
  MakeNNStarToNN<channelType, N1675pPC, N1675nPC>::Make(this);
  MakeNNStarToNN<channelType, N1680pPC, N1680nPC>::Make(this);
  MakeNNStarToNN<channelType, N1700pPC, N1700nPC>::Make(this);
  MakeNNStarToNN<channelType, N1710pPC, N1710nPC>::Make(this);
  MakeNNStarToNN<channelType, N1720pPC, N1720nPC>::Make(this);
  MakeNNStarToNN<channelType, N1900pPC, N1900nPC>::Make(this);
  MakeNNStarToNN<channelType, N1990pPC, N1990nPC>::Make(this);
  MakeNNStarToNN<channelType, N2090pPC, N2090nPC>::Make(this);
  MakeNNStarToNN<channelType, N2190pPC, N2190nPC>::Make(this);
  MakeNNStarToNN<channelType, N2220pPC, N2220nPC>::Make(this);
  MakeNNStarToNN<channelType, N2250pPC, N2250nPC>::Make(this);
}

