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
// $Id: G4CollisionNNToNNstar.cc,v 1.2 2003-12-12 13:34:16 hpw Exp $ //

#include "globals.hh"
#include "G4CollisionNNToNNstar.hh"
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
#include "G4CollisionNNToNNstar.hh"
#include "G4ConcreteNNToNNStar.hh"
#include "G4HadParticleCodes.hh"
#include "G4HadParticleCodes.hh"
#include "G4Pair.hh"

typedef G4ConcreteNNToNNStar channelType;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N1400nPC)  theC1;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N1400pPC)  theC2;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N1400pPC)  theC3;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N1400nPC)  theC4;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N1520nPC)  theC5;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N1520pPC)  theC6;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N1520pPC)  theC7;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N1520nPC)  theC8;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N1535nPC)  theC9;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N1535pPC)  theC10;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N1535pPC)  theC11;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N1535nPC)  theC12;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N1650nPC)  theC13;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N1650pPC)  theC14;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N1650pPC)  theC15;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N1650nPC)  theC16;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N1675nPC)  theC17;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N1675pPC)  theC18;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N1675pPC)  theC19;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N1675nPC)  theC20;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N1680nPC)  theC21;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N1680pPC)  theC22;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N1680pPC)  theC23;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N1680nPC)  theC24;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N1700nPC)  theC25;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N1700pPC)  theC26;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N1700pPC)  theC27;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N1700nPC)  theC28;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N1710nPC)  theC29;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N1710pPC)  theC30;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N1710pPC)  theC31;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N1710nPC)  theC32;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N1720nPC)  theC33;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N1720pPC)  theC34;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N1720pPC)  theC35;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N1720nPC)  theC36;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N1900nPC)  theC37;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N1900pPC)  theC38;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N1900pPC)  theC39;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N1900nPC)  theC40;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N1990nPC)  theC41;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N1990pPC)  theC42;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N1990pPC)  theC43;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N1990nPC)  theC44;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N2090nPC)  theC45;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N2090pPC)  theC46;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N2090pPC)  theC47;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N2090nPC)  theC48;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N2190nPC)  theC49;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N2190pPC)  theC50;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N2190pPC)  theC51;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N2190nPC)  theC52;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N2220nPC)  theC53;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N2220pPC)  theC54;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N2220pPC)  theC55;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N2220nPC)  theC56;

typedef INT4(channelType, NeutronPC, NeutronPC, NeutronPC, N2250nPC)  theC57;
typedef INT4(channelType, ProtonPC,  ProtonPC,  ProtonPC,  N2250pPC)  theC58;
typedef INT4(channelType, NeutronPC, ProtonPC,  NeutronPC, N2250pPC)  theC59;
typedef INT4(channelType, NeutronPC, ProtonPC,  ProtonPC,  N2250nPC)  theC60;


typedef GROUP60(theC1, theC2, theC3, theC4, theC5, theC6, theC7, theC8, theC9, theC10,
              theC11, theC12, theC13, theC14, theC15, theC16, theC17, theC18, theC19, theC20,
              theC21, theC22, theC23, theC24, theC25, theC26, theC27, theC28, theC29, theC30,
              theC31, theC32, theC33, theC34, theC35, theC36, theC37, theC38, theC39, theC40,
              theC41, theC42, theC43, theC44, theC45, theC46, theC47, theC48, theC49, theC50,
              theC51, theC52, theC53, theC54, theC55, theC56, theC57, theC58, theC59, theC60) theChannels;

G4CollisionNNToNNstar::G4CollisionNNToNNstar()
{ 
  G4ForEach<theChannels, Resolve>::Apply(this);
}
