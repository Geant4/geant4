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
#include "G4Neutron.hh"
#include "G4NeutronHPElasticData.hh"


main()
{
   G4Element* theElement = new G4Element("carbon", "C", 6, 12.0*g/mole);
//   G4Element* theElement = new G4Element("Plutonium", "Pu", 1);
//   G4Isotope * theIso = new G4Isotope("Pu140", 94., 240, 240.0*g/mole);
//   theElement->AddIsotope(theIso, 100.*perCent);

   G4NeutronHPElasticData anElasticDataSet; 
   G4ParticleDefinition* theParticleDefinition = G4Neutron::NeutronDefinition();

   G4double temp;
   G4cin >> temp;
   temp *= kelvin;
//   G4double ekin = 0.5E+04*eV;
   G4double ekin = 1.E-5*eV;
//   G4double ekin = 1.E-3*eV;
   G4DynamicParticle* theDynamicParticle;
   G4int counter = -1;
   G4int points = -1; 
   G4int hpw=0;
   while (++points<20)
//   while (++points<500)
   {
//     ekin +=0.002E+04*eV;
     ekin *=3.333333333; 
//     if(ekin<0.8*eV) ekin *=1.2;
//     else if(ekin<12*eV) ekin*=1.05;
//     else ekin*=1.005;
     theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(1.,0.,0.), ekin);
     while(++counter<50)
//     while(++counter<1)
     {
       hpw++;
       if(hpw == 1*(hpw/1)) G4cerr << "point number "<<hpw<<" being processed"<<endl;
       cout << ekin/MeV 
            << " " 
            << anElasticDataSet.GetCrossSection(theDynamicParticle, theElement, temp)/millibarn
            << G4endl;
     }
   delete theDynamicParticle;
   counter = -1;
   }
}
