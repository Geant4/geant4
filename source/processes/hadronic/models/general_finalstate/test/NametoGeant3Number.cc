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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: NametoGeant3Number.cc,v 1.1 2003-10-08 12:04:17 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 G4int NameToGeant3Number( G4ParticleDefinition* p )
  {
    G4String name = p->GetParticleName();
    G4int n = 0;
    if( name == "gamma" )
      n = 0;
    else if( name == "positron" )
      n = 2;
    else if( name == "electron" )
      n = 3;
    else if( name == "neutrino" )
      n = 4;
    else if( name == "mu+" )
      n = 5;
    else if( name == "mu-" )
      n = 6;
    else if( name == "pi0" )
      n = 7;
    else if( name == "pi+" )
      n = 8;
    else if( name == "pi-" )
      n = 9;
    else if( name == "kaon0L" )
      n = 10;
    else if( name == "kaon+" )
      n = 11;
    else if( name == "kaon-" )
      n = 12;
    else if( name == "neutron" )
      n = 13;
    else if( name == "proton" )
      n = 14;
    else if( name == "anti_proton" )
      n = 15;
    else if( name == "kaon0S" )
      n = 16;
    else if( name == "lambda" )
      n = 18;
    else if( name == "sigma+" )
      n = 19;
    else if( name == "sigma0" )
      n = 20;
    else if( name == "sigma-" )
      n = 21;
    else if( name == "xi0" )
      n = 22;
    else if( name == "xi-" )
      n = 23;
    else if( name == "omega-" )
      n = 24;
    else if( name == "anti_neutron" )
      n = 25;
    else if( name == "anti_lambda" )
      n = 26;
    else if( name == "anti_sigma-" )
      n = 27;
    else if( name == "anti_sigma0" )
      n = 28;
    else if( name == "anti_sigma+" )
      n = 29;
    else if( name == "anti_xi0" )
      n = 30;
    else if( name == "anti_xi-" )
      n = 31;
    else if( name == "anti_omega-" )
      n = 32;
    else if( name == "deuteron" )
      n = 45;
    else if( name == "triton" )
      n = 46;
    else if( name == "alpha" )
      n = 47;
    return n;
  }
