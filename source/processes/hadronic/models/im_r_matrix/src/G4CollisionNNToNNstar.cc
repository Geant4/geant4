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
// $Id: G4CollisionNNToNNstar.cc,v 1.3 2003-12-12 15:38:22 hpw Exp $ //

#include "globals.hh"
#include "G4CollisionNNToNNstar.hh"
#include "G4ConcreteNNToNNStar.hh"
#include "G4HadParticleCodes.hh"
#include "G4Pair.hh"

typedef G4ConcreteNNToNNStar channelType;

G4CollisionNNToNNstar::G4CollisionNNToNNstar()
{ 
  MakeNNToNNStar<N1400pPC, N1400nPC, channelType>(this);
  MakeNNToNNStar<N1520pPC, N1520nPC, channelType>(this);
  MakeNNToNNStar<N1535pPC, N1535nPC, channelType>(this);
  MakeNNToNNStar<N1650pPC, N1650nPC, channelType>(this);
  MakeNNToNNStar<N1675pPC, N1675nPC, channelType>(this);
  MakeNNToNNStar<N1680pPC, N1680nPC, channelType>(this);
  MakeNNToNNStar<N1700pPC, N1700nPC, channelType>(this);
  MakeNNToNNStar<N1710pPC, N1710nPC, channelType>(this);
  MakeNNToNNStar<N1720pPC, N1720nPC, channelType>(this);
  MakeNNToNNStar<N1900pPC, N1900nPC, channelType>(this);
  MakeNNToNNStar<N2090pPC, N2090nPC, channelType>(this);
  MakeNNToNNStar<N2190pPC, N2190nPC, channelType>(this);
  MakeNNToNNStar<N2220pPC, N2220nPC, channelType>(this);
  MakeNNToNNStar<N2250pPC, N2250nPC, channelType>(this);
}
