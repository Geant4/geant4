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
// $Id: G4AntiNucleiConstructor.cc,v 1.1 2010-06-11 05:50:20 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      first implementaion                11  June 2010 H.Kurashige
//

#include "G4ParticleTable.hh"
#include "G4AntiNucleus.hh"
#include "G4AntiNucleiConstructor.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"


G4AntiNucleiConstructor::G4AntiNucleiConstructor()
{
}

G4AntiNucleiConstructor::~G4AntiNucleiConstructor()
{
}

void G4AntiNucleiConstructor::Construct()
{   
  // create particle
  //
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived     
  G4ParticleDefinition* p;
  
  // Deuteron
  p = new G4AntiNucleus(
      "anti_deuteron",    1.875613*GeV,       0.0*MeV,  -1.0*eplus,
                    2,              +1,             0,
                    0,               0,             0,
       "anti_nucleus",               0,            -2, -1000010020,
                 true,            -1.0,          NULL   
      );
  // Magnetic Moment
  G4double mN = eplus*hbar_Planck/2./(proton_mass_c2 /c_squared);
  p->SetPDGMagneticMoment( 0.857438230 * mN);

  //particle registeration
  p->SetAntiPDGEncoding(1000010020);


  // Triton
  p = new G4AntiNucleus(
        "anti_triton",    2.808921*GeV,       0.0*MeV,  -1.0*eplus,
                    1,              +1,             0,
                    0,               0,             0,
       "anti_nucleus",               0,            -3, -1000010030,
                 true,            -1.0,          NULL
     );
  
  // Magnetic Moment
  mN = eplus*hbar_Planck/2./(proton_mass_c2 /c_squared);
  p->SetPDGMagneticMoment( 2.97896248 * mN);

  //particle registeration
  p->SetAntiPDGEncoding(1000010030);


  // He3
  p = new G4AntiNucleus(
           "anti_He3",    2.808391*GeV,       0.0*MeV,  -2.0*eplus,
                    1,              +1,             0,
                    0,               0,             0,
       "anti_nucleus",               0,            -3, -1000020030,
                 true,            -1.0,          NULL
	   );
  
  // Magnetic Moment
  mN = eplus*hbar_Planck/2./(proton_mass_c2 /c_squared);
  p->SetPDGMagneticMoment( 2.97896248 * mN);

  //particle registeration
  p->SetAntiPDGEncoding(1000020030);

  // He4 (Alpha)
  p = new G4AntiNucleus(
	 "anti_alpha",    3.727379*GeV,       0.0*MeV,  -2.0*eplus,
                    0,              +1,             0,
                    0,               0,             0,
       "anti_nucleus",               0,            -4,  -1000020040,
                 true,            -1.0,          NULL
	   );
  
  // Magnetic Moment
  mN = eplus*hbar_Planck/2./(proton_mass_c2 /c_squared);
  p->SetPDGMagneticMoment( 2.97896248 * mN);

  //particle registeration
  p->SetAntiPDGEncoding(1000020040);

}
