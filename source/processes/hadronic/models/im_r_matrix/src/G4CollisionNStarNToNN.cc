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
// $Id: G4CollisionNStarNToNN.cc,v 1.2.2.1 2004/03/24 13:18:44 hpw Exp $ //

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

