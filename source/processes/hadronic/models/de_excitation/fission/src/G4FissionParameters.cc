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
// $Id: G4FissionParameters.cc 89550 2015-04-17 08:38:15Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
//J. M. Quesada (May 2009): sigma_sym (SigmaS) tuned for spallation data.
//J. M. Quesada (30.10.09): retuning for IAEA spallation data

#include "G4FissionParameters.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exp.hh"

G4FissionParameters::G4FissionParameters()
  : A1(134),A2(141),A3((A1 + A2)*0.5),As(0.0),Sigma1(0.0),Sigma2(0.0),
    SigmaS(0.),w(0.0)
{}

G4FissionParameters::~G4FissionParameters()
{}

void G4FissionParameters::DefineParameters(G4int A, G4int Z, G4double ExEnergy,
					   G4double FissionBarrier)
{
  // to avoid usage of units
  G4double U = ExEnergy/CLHEP::MeV; 

  As = A*0.5;
    
  if (A <= 235) { Sigma2 = 5.6; }  
  else          { Sigma2 = 5.6 + 0.096*(A-235); }

  Sigma1 = 0.5*Sigma2;     
    
  //JMQ 310509 
  //    if (SigmaS > 20.0) SigmaS = 20.0;
  //   SigmaS*=1.3;
  //JMQ 301009: retuning (after CEM transition prob.have been chosen as default)
  SigmaS = 0.8*G4Exp(0.00553*U + 2.1386); 
    
  G4double x1 = (A1-As)/Sigma1;
  G4double x2 = (A2-As)/Sigma2;
  G4double FasymAsym = 2*G4Exp(-0.5*x2*x2) + G4Exp(-0.5*x1*x1);
    
  G4double x3 = (As-A3)/SigmaS;
  G4double FsymA1A2 = G4Exp(-0.5*x3*x3);
    
  G4double wa = 0.0;
  w = 0.0;
  if (Z >= 90) {  
    if (U <= 16.25) { wa = G4Exp(0.5385*U-9.9564); }  
    else            { wa = G4Exp(0.09197*U-2.7003); }            
  } else if (Z == 89) {  
    wa = G4Exp(0.09197*U-1.0808);
  } else if (Z >= 82) {  
    G4double X = std::max(0.0, FissionBarrier/CLHEP::MeV - 7.5);
    wa = G4Exp(0.09197*(U-X) - 1.0808);
  } else {               // Z < 82
    w = 1001.0;
  }
    
  if (w == 0.0) {
    G4double w1 = std::max(1.03*wa - FasymAsym, 0.0001);
    G4double w2 = std::max(1.0 - FsymA1A2*wa,   0.0001);
    
    w = w1/w2;
      
    if (82 <= Z && Z < 89 && A < 227) { w *= G4Exp(0.3*(227-A)); }
  }
}


