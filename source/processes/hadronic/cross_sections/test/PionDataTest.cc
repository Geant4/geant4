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
   G4cout << "Thank you !"<<G4endl<<G4endl;
   G4Element* theElement = new G4Element("El    ", "any",  z,  a*g/mole);
   
   G4cout << "Please choose the particle type: 1=Pi+, 2=Pi-"<<G4endl;
   G4int ptype;
   G4cin >> ptype;
   G4cout << "Thank you !"<<G4endl<<G4endl;
   
   G4ParticleDefinition* theParticleDefinition = G4PionPlus::PionPlusDefinition();
   if(2==ptype) theParticleDefinition = G4PionMinus::PionMinusDefinition();

   G4double ekin = 1*MeV;
   G4DynamicParticle* theDynamicParticle;
   G4int count = 0;
   while(ekin<10*GeV)
   {
     ekin *= 1.2;
     theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                                 G4ParticleMomentum(1.,0.,0.), ekin);
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
