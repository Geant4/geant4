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
// $Id: G4CollisionNNToDeltaNstar.cc,v 1.3.2.1 2004/03/24 13:18:33 hpw Exp $ //

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

