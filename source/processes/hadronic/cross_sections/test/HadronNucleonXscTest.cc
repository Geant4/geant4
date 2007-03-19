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
// $Id: HadronNucleonXscTest.cc,v 1.1 2007-03-19 10:45:14 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Test of G4CrossSectionDataStore and G4CrossSectionDataSet classes
//
// History:
//
// 29.05.06 V.Grichine: NIST elements/materials, write in file
// F.W. Jones, TRIUMF, 22-JAN-98
//                     19-MAY-98
//



#include "G4ios.hh"
#include <fstream>
#include <cmath>
#include "globals.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"

#include <iomanip>



#include "G4PionPlus.hh"
#include "G4Proton.hh"

#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"



#include "G4HadronNucleonXsc.hh"




#ifdef G4_SOLVE_TEMPLATES
#include "g4templates.hh"
#endif


int main()
{
  G4int choice;
// Particle definition

  G4cout << " 1 proton" << G4endl;
  G4cout << " 2 pion+" << G4endl;
  G4cout << " 3 neutron" << G4endl;
  G4cout << " 4 kaon+" << G4endl;
  G4cout << " 5 kaon0short" << G4endl;
  G4cout << " 6 pion-" << G4endl;
  //  G4cin >> choice;
  choice = 3;

  G4ParticleDefinition* theParticleDefinition;

  switch (choice) 
  {
    case 1:

      theParticleDefinition = G4Proton::ProtonDefinition();
      break;

    case 2:

      theParticleDefinition = G4PionPlus::PionPlusDefinition();
      break;

    case 3:

      theParticleDefinition = G4Neutron::NeutronDefinition();
      break;

    case 4:

      theParticleDefinition = G4KaonPlus::KaonPlusDefinition();
      break;

    case 5:

      theParticleDefinition = G4KaonZeroShort::KaonZeroShortDefinition();
      break;

    case 6:

      theParticleDefinition = G4PionMinus::PionMinusDefinition();
      break;
  }
  // Nucleon definition

  G4cout << " 1 proton" << G4endl;
  G4cout << " 2 neutron" << G4endl;
  choice = 1;

  G4ParticleDefinition* theNucleon;

  switch (choice) 
  {
    case 1:

      theNucleon = G4Proton::ProtonDefinition();
      break;

    case 3:

      theNucleon = G4Neutron::NeutronDefinition();
      break;

  }
  // Parametrisation definition

  G4cout << " 1 EL" << G4endl;
  G4cout << " 2 PDG" << G4endl;
  G4cout << " 3 NS" << G4endl;
  G4cout << " 4 VU" << G4endl;
  choice = 1;


  G4int i, iMax;
  G4double kinEnergy;
 
  G4DynamicParticle* theDynamicParticle;
  G4HadronNucleonXsc* hnXsc = new G4HadronNucleonXsc();

  G4double sig = 0;  

  std::ofstream writef("hnxsc.dat", std::ios::out ) ;
  writef.setf( std::ios::scientific, std::ios::floatfield );

  kinEnergy = 0.01*GeV;
  iMax = 90;
    
  writef <<iMax<< G4endl; 

  for(i = 0; i < iMax; i++)
  {
   
    theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(1.,0.,0.), 
                                              kinEnergy);


  switch (choice) 
  {
    case 1:
      sig = hnXsc->GetHadronNucleonXscEL(theDynamicParticle,theNucleon);
    break;

    case 2:
      sig = hnXsc->GetHadronNucleonXscPDG(theDynamicParticle,theNucleon);
    break;

    case 3:
      sig = hnXsc->GetHadronNucleonXscNS(theDynamicParticle,theNucleon);
    break;

    case 4:
      sig = hnXsc->GetHadronNucleonXscVU(theDynamicParticle,theNucleon);
    break;
  }
      

  G4cout << kinEnergy/GeV << " GeV, \t"<< " UV xsc = " <<" \t"<< sig/millibarn << " mb" << G4endl;


    writef << kinEnergy/GeV <<"\t"<< sig/millibarn << G4endl;

    kinEnergy *= 1.138;
    delete theDynamicParticle;
  }


  G4cout<<"energy in GeV"<<"\t"<<"cross-section in millibarn"<<G4endl;
  G4cout << " cross section for " << 
            theParticleDefinition->GetParticleName() <<
           " on " << theNucleon->GetParticleName() << G4endl;
  


  return 1;
} // end of main
