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
// $Id: G4FissionParameters.cc,v 1.8 2010-11-17 20:22:46 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
//J. M. Quesada (May 2009): sigma_sym (SigmaS) tuned for spallation data.
//J. M. Quesada (30.10.09): retuning for IAEA spallation data

#include "G4FissionParameters.hh"
#include "G4HadronicException.hh"


const G4double G4FissionParameters::A1 = 134.0;
const G4double G4FissionParameters::A2 = 141.0;

G4FissionParameters::G4FissionParameters(G4int A, G4int Z, G4double ExEnergy,
					 G4double FissionBarrier)
{
  G4double U = ExEnergy; 
  
  As = A/2.0;
    
  if (A <= 235) Sigma2 = 5.6;  // MeV
  else Sigma2 = 5.6 + 0.096*(A-235); // MeV
    
  Sigma1 = 0.5*Sigma2; // MeV
    
  SigmaS = std::exp(0.00553*U/MeV + 2.1386); // MeV
    
  //JMQ 310509 
  //    if (SigmaS > 20.0) SigmaS = 20.0;
  //   SigmaS*=1.3;
  //JMQ 301009: retuning (after CEM transition prob.have been chosen as default)
  SigmaS*=0.8;
  //
    
  G4double FasymAsym = 2.0*std::exp(-((A2-As)*(A2-As))/(2.0*Sigma2*Sigma2)) + 
    std::exp(-((A1-As)*(A1-As))/(2.0*Sigma1*Sigma1));
    
  G4double FsymA1A2 = std::exp(-((As-(A1+A2)/2.0)*(As-(A1+A2)/2.0))
			       /(2.0*SigmaS*SigmaS));
    
    
  G4double wa = 0.0;
  w = 0.0;
  if (Z >= 90) {         // Z >= 90
    if (U <= 16.25) wa = std::exp(0.5385*U/MeV-9.9564);  // U <= 16.25 MeV
    else wa = std::exp(0.09197*U/MeV-2.7003);            // U  > 16.25 MeV
  } else if (Z == 89) {  // Z == 89
    wa = std::exp(0.09197*U-1.0808);
  } else if (Z >= 82) {  //  82 <= Z <= 88
    G4double X = FissionBarrier - 7.5*MeV;
    if (X < 0.0) X = 0.0;
    wa = std::exp(0.09197*(U-X)/MeV-1.0808);
  } else {               // Z < 82
    w = 1001.0;
  }
    
  if (w == 0.0) {
    G4double w1 = std::max(1.03*wa - FasymAsym, 0.0001);
    G4double w2 = std::max(1.0 - FsymA1A2*wa,   0.0001);
    
    w = w1/w2;
      
    if (82 <= Z && Z < 89 && A < 227)  w *= std::exp(0.3*(227-A));
  }
    
}

G4FissionParameters::~G4FissionParameters()
{}

