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
// $Id: G4CollisionNNToDeltaNstar.cc,v 1.4 2010-03-12 15:45:18 gunter Exp $ //

#include "globals.hh"
#include "G4CollisionNNToDeltaNstar.hh"
#include "G4ConcreteNNToDeltaNstar.hh"
#include "G4HadParticleCodes.hh"
#include "G4Pair.hh"

typedef G4ConcreteNNToDeltaNstar channelType;

G4CollisionNNToDeltaNstar::G4CollisionNNToDeltaNstar()
{ 
  MakeNNToDeltaNstar<N1400pPC, channelType, N1400nPC>::Make(this);
  MakeNNToDeltaNstar<N1520pPC, channelType, N1520nPC>::Make(this);
  MakeNNToDeltaNstar<N1535pPC, channelType, N1535nPC>::Make(this);
  MakeNNToDeltaNstar<N1650pPC, channelType, N1650nPC>::Make(this);
  MakeNNToDeltaNstar<N1675pPC, channelType, N1675nPC>::Make(this);
  MakeNNToDeltaNstar<N1680pPC, channelType, N1680nPC>::Make(this);
  MakeNNToDeltaNstar<N1700pPC, channelType, N1700nPC>::Make(this);
  MakeNNToDeltaNstar<N1710pPC, channelType, N1710nPC>::Make(this);
  MakeNNToDeltaNstar<N1720pPC, channelType, N1720nPC>::Make(this);
  MakeNNToDeltaNstar<N1900pPC, channelType, N1900nPC>::Make(this);
  MakeNNToDeltaNstar<N1990pPC, channelType, N1990nPC>::Make(this);
  MakeNNToDeltaNstar<N2090pPC, channelType, N2090nPC>::Make(this);
  MakeNNToDeltaNstar<N2190pPC, channelType, N2190nPC>::Make(this);
  MakeNNToDeltaNstar<N2220pPC, channelType, N2220nPC>::Make(this);
  MakeNNToDeltaNstar<N2250pPC, channelType, N2250nPC>::Make(this);
}

