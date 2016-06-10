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
// $Id: G4N14GEMProbability.cc 87017 2014-11-21 16:26:26Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4N14GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4N14GEMProbability::G4N14GEMProbability() :
  G4GEMProbability(14,7,1.0) // A,Z,Spin
{

  ExcitEnergies.push_back(2312.798*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(60e-15*s);
	  
  ExcitEnergies.push_back(3948.1*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(3.1e-15*s);
		
  ExcitEnergies.push_back(4915.1*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(5.3e-15*s);
		
  ExcitEnergies.push_back(5105.89*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(8.6e-12*s);
	  
  ExcitEnergies.push_back(5691.44*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(6.9e-15*s);
		
  ExcitEnergies.push_back(5834.25*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(12.5e-12*s);
		
  ExcitEnergies.push_back(6203.5*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(86e-15*s);
		
  ExcitEnergies.push_back(6446.17*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(435e-15*s);
		
  ExcitEnergies.push_back(7029.12*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(3.7e-15*s);
		
  ExcitEnergies.push_back(7966.9*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(2.5E-3*keV));
		
  ExcitEnergies.push_back(8490.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(27E-6*keV));
		
  ExcitEnergies.push_back(8618*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(7.0*keV));
		
  ExcitEnergies.push_back(8776*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(460*keV));
		
  ExcitEnergies.push_back(8907.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(19.7*keV));
		
  ExcitEnergies.push_back(8964.0 *keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(7E-6*keV));
	  
  ExcitEnergies.push_back(8979*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(8*keV));
		  
  ExcitEnergies.push_back(9129*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(1*keV));
		  
  ExcitEnergies.push_back(9172.25*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(0.074*keV));
		
  ExcitEnergies.push_back(9386.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(15.6 *keV));
		
  ExcitEnergies.push_back(9509*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(41*keV));
		  
  ExcitEnergies.push_back(9703 *keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(15*keV));
		 
  ExcitEnergies.push_back(10063*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(10*keV));
		
  ExcitEnergies.push_back(10101*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(5*keV));
		
  ExcitEnergies.push_back(10226*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(80*keV));
		
  ExcitEnergies.push_back(10432*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(33*keV));
		
  ExcitEnergies.push_back(10540*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(140*keV));
		
  ExcitEnergies.push_back(10812*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(0.39*eV));
		
  ExcitEnergies.push_back(11050*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(2*keV));
		
  ExcitEnergies.push_back(11070*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(100*keV));
		
  ExcitEnergies.push_back(11240*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(20*keV));
		
  ExcitEnergies.push_back(11290*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(180*keV));
		
  ExcitEnergies.push_back(11357*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(30*keV));
		
  ExcitEnergies.push_back(11513.6*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(7.0*keV));
	  
  ExcitEnergies.push_back(11680*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(150*keV));
		
  ExcitEnergies.push_back(11741*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(40*keV));
		
  ExcitEnergies.push_back(11761*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(78*keV));
		
  ExcitEnergies.push_back(11807*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(119*keV));
		
  ExcitEnergies.push_back(11874*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(101*keV));
		
  ExcitEnergies.push_back(12200*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(300*keV));
		
  ExcitEnergies.push_back(12408*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(37*keV));
		
  ExcitEnergies.push_back(12418*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(41*keV));
		
  ExcitEnergies.push_back(12594*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(48*keV));
		
  ExcitEnergies.push_back(12688*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(22*keV));
		
  ExcitEnergies.push_back(12792*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(18*keV));
		
  ExcitEnergies.push_back(12819*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(8*keV));
		
  ExcitEnergies.push_back(12923*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(25*keV));
		
  ExcitEnergies.push_back(13166*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(15*keV));
		
  ExcitEnergies.push_back(13192*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(65*keV));
		
  ExcitEnergies.push_back(13243*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(92*keV));
		
  ExcitEnergies.push_back(13300*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(1000*keV));
		
  ExcitEnergies.push_back(13656*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(90*keV));
		
  ExcitEnergies.push_back(13714*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(105*keV));
		
  ExcitEnergies.push_back(13710*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(180*keV));
		
  ExcitEnergies.push_back(14250*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(420*keV));
		
  ExcitEnergies.push_back(14660*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(100*keV));
		
  ExcitEnergies.push_back(16800*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(300*keV));
		
  ExcitEnergies.push_back(16910*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(100*keV));
		
  ExcitEnergies.push_back(17170*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(300*keV));
		
  ExcitEnergies.push_back(18100*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(300*keV));
		
  ExcitEnergies.push_back(18100*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(600*keV));
		
  ExcitEnergies.push_back(18200*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(400*keV));
		
  ExcitEnergies.push_back(18400*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(300*keV));
		
  ExcitEnergies.push_back(18500*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(60 *keV));
		
  ExcitEnergies.push_back(18800*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(400*keV));
		
  ExcitEnergies.push_back(20100*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(500*keV));
		
  ExcitEnergies.push_back(20800*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(600*keV));
		
  ExcitEnergies.push_back(20800*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(500*keV));
		
  ExcitEnergies.push_back(21300*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(1000*keV));
		
  ExcitEnergies.push_back(21500*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(500*keV));
		
  ExcitEnergies.push_back(21700*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(200*keV));
		
  ExcitEnergies.push_back(23000*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(3000*keV));
		
  ExcitEnergies.push_back(23.3E3*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(500*keV));
}

G4N14GEMProbability::~G4N14GEMProbability() 
{}
