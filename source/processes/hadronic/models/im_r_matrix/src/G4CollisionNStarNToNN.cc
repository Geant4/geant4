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
// $Id: G4CollisionNStarNToNN.cc,v 1.2 2003-12-15 10:12:48 hpw Exp $ //

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
  MakeNNStarToNN<channelType, N1400pPC, N1400nPC>(this);
  MakeNNStarToNN<channelType, N1520pPC, N1520nPC>(this);
  MakeNNStarToNN<channelType, N1535pPC, N1535nPC>(this);
  MakeNNStarToNN<channelType, N1650pPC, N1650nPC>(this);
  MakeNNStarToNN<channelType, N1675pPC, N1675nPC>(this);
  MakeNNStarToNN<channelType, N1680pPC, N1680nPC>(this);
  MakeNNStarToNN<channelType, N1700pPC, N1700nPC>(this);
  MakeNNStarToNN<channelType, N1710pPC, N1710nPC>(this);
  MakeNNStarToNN<channelType, N1720pPC, N1720nPC>(this);
  MakeNNStarToNN<channelType, N1900pPC, N1900nPC>(this);
  MakeNNStarToNN<channelType, N1990pPC, N1990nPC>(this);
  MakeNNStarToNN<channelType, N2090pPC, N2090nPC>(this);
  MakeNNStarToNN<channelType, N2190pPC, N2190nPC>(this);
  MakeNNStarToNN<channelType, N2220pPC, N2220nPC>(this);
  MakeNNStarToNN<channelType, N2250pPC, N2250nPC>(this);
}

