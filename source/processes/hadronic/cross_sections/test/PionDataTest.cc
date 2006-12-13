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
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4HadronCrossSections.hh"


main()
{
   G4PiNuclearCrossSection aSpecialDataSet;
   G4HadronCrossSections aGeneralDataSet;
   
   G4cout << "Please enter element z, a"<<G4endl;
   G4int z;
   G4double a;
   G4cin >> z >> a;
   G4cout << "Thank you !"<<G4endl;
   G4Element* theElement = new G4Element("El    ", "any",  z,  a*g/mole);
   
   G4cout << "Please choose the particle type: 1=Pi+, 2=Pi-"<<G4endl;
   G4int ptype;
   G4cin >> ptype;
   G4cout << "Thank you !"<<G4endl;
   
   G4ParticleDefinition* theParticleDefinition = G4PionPlus::PionPlusDefinition();
   if(2==ptype) theParticleDefinition = G4PionMinus::PionMinusDefinition();

   G4double ekin = 1*MeV;
   G4DynamicParticle* theDynamicParticle;
   G4int count = 0;
   while(ekin<10*GeV)
   {
     ekin *= 1.02;
     theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                                 G4ParticleMomentum(1.,0.,0.), ekin);
if( ! aSpecialDataSet.IsApplicable(theDynamicParticle, theElement))
{
  std::cout << "No way we use this here "<<std::endl;
  abort();
}

//     if(aDataSet.IsApplicable(theDynamicParticle, theElement))
//     {
       G4cout << ekin/GeV  << " " 
            << aGeneralDataSet.GetInelasticCrossSection(theDynamicParticle, theElement)/millibarn<<" "
            << aSpecialDataSet.GetCrossSection(theDynamicParticle, theElement, 273*kelvin)/millibarn
            << G4endl;
//     }
     delete theDynamicParticle;
   }
}
