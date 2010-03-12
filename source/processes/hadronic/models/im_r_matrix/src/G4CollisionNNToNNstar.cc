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
// $Id: G4CollisionNNToNNstar.cc,v 1.5 2010-03-12 15:45:18 gunter Exp $ //

#include "globals.hh"
#include "G4CollisionNNToNNstar.hh"
#include "G4ConcreteNNToNNStar.hh"
#include "G4HadParticleCodes.hh"
#include "G4Pair.hh"

typedef G4ConcreteNNToNNStar channelType;

G4CollisionNNToNNstar::G4CollisionNNToNNstar()
{ 
  MakeNNToNNStar<N1400pPC, N1400nPC, channelType>::Make(this);
  MakeNNToNNStar<N1520pPC, N1520nPC, channelType>::Make(this);
  MakeNNToNNStar<N1535pPC, N1535nPC, channelType>::Make(this);
  MakeNNToNNStar<N1650pPC, N1650nPC, channelType>::Make(this);
  MakeNNToNNStar<N1675pPC, N1675nPC, channelType>::Make(this);
  MakeNNToNNStar<N1680pPC, N1680nPC, channelType>::Make(this);
  MakeNNToNNStar<N1700pPC, N1700nPC, channelType>::Make(this);
  MakeNNToNNStar<N1710pPC, N1710nPC, channelType>::Make(this);
  MakeNNToNNStar<N1720pPC, N1720nPC, channelType>::Make(this);
  MakeNNToNNStar<N1900pPC, N1900nPC, channelType>::Make(this);
  MakeNNToNNStar<N1990pPC, N1990nPC, channelType>::Make(this);
  MakeNNToNNStar<N2090pPC, N2090nPC, channelType>::Make(this);
  MakeNNToNNStar<N2190pPC, N2190nPC, channelType>::Make(this);
  MakeNNToNNStar<N2220pPC, N2220nPC, channelType>::Make(this);
  MakeNNToNNStar<N2250pPC, N2250nPC, channelType>::Make(this);
}
