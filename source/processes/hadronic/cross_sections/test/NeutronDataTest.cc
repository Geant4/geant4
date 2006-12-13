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
#include "G4Neutron.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4Proton.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4HadronCrossSections.hh"
#include "G4NeutronHPCaptureData.hh"


main()
{
   G4NeutronInelasticCrossSection aNDataSet;
   G4ProtonInelasticCrossSection aPDataSet;
   G4HadronCrossSections aGeneralDataSet;
   G4NeutronHPCaptureData aLowNDataSet;
//   G4Element* theElement = new G4Element("copper", "Cu", 29, 63.54*g/mole);
//   G4Element* theElement = new G4Element("copper", "Al", 13, 27.0*g/mole);
//   G4Element* theElement = new G4Element("be    ", "Be",  4,  9.0*g/mole);
   G4Element* theElement = new G4Element("H    ", "H",  1,  1.0*g/mole);
   G4ParticleDefinition* theParticleDefinition = G4Neutron::NeutronDefinition();

   G4double ekin = 0.0001*eV;
   G4DynamicParticle* theDynamicParticle;
   G4int count = 0;
   while(ekin<10*keV)
   {
     ekin *= 1.2;
     theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                                 G4ParticleMomentum(1.,0.,0.), ekin);
//     if(aDataSet.IsApplicable(theDynamicParticle, theElement))
     {
       cout << ekin/eV 
            << " " 
            << aLowNDataSet.GetCrossSection(theDynamicParticle, theElement, 273*kelvin)/millibarn
            << G4endl;
     }
     delete theDynamicParticle;
   }
}
