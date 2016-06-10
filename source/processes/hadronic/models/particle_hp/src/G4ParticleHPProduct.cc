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
// particle_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 080718 As for secondary photons, if its mean value has a value of integer, 
//        then a sampling of multiplicity that based on Poisson Distribution 
//        is not carried out and the mean is used as a multiplicity.
//        modified by T. Koi.
// 080721 Using ClearHistories() methodl for limiting the sum of secondary energies
//        modified by T. Koi.
// 080901 bug fix of too many secnodaries production in nd reactinos by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPProduct.hh" 
#include "G4Poisson.hh"
#include "G4Proton.hh"

G4int G4ParticleHPProduct::GetMultiplicity(G4double anEnergy )
{
  //if(theDist == 0) { return 0; }
  //151120 TK Modified for solving reproducibility problem 
  if ( theDist == 0 ) { 
     fCache.Get().theCurrentMultiplicity = 0;
     return 0; 
  }

  G4double mean = theYield.GetY(anEnergy);
  //g  G4cout << "G4ParticleHPProduct MEAN NUMBER OF PARTICLES " << mean << " for " << theMass << G4endl;
  //if( mean <= 0. ) return 0;
  //151120 TK Modified for solving reproducibility problem 
  //This is also a real fix
  if ( mean <= 0. )  {
     fCache.Get().theCurrentMultiplicity = 0;
     return 0; 
  }

  G4int multi;
  multi = G4int(mean+0.0001);
  //if(theMassCode==0) multi = G4Poisson(mean); // @@@@gammas. please X-check this
  //080718
#ifdef PHP_AS_HP
  if ( theMassCode == 0 ) // DELETE THIS: IT MUST BE DONE FOR ALL PARTICLES
#endif
  { 
     if ( G4int ( mean ) == mean )
     {
        multi = (G4int) mean;
     }
     else
     {
#ifdef PHP_AS_HP
	 multi = G4Poisson ( mean ); 
#else 
       if( theMultiplicityMethod == G4HPMultiPoisson ) {
	 multi = G4Poisson ( mean ); 	 
	 if( getenv("G4PHPTEST") )  G4cout << " MULTIPLICITY MULTIPLIED " << multi << " " << theMassCode << G4endl;
       } else { // if( theMultiplicityMethod == G4HPMultiBetweenInts ) {
	 G4double radnf = CLHEP::RandFlat::shoot();
	 G4int imulti = G4int(mean);
	 multi = imulti + G4int(radnf < mean-imulti);
	 //	 G4cout << theMass << " multi " << multi << " mean " << mean 
	 //		<< " radnf " << radnf << " mean-imulti " << mean-imulti << G4endl;
       }
#endif
	//       multi = int(mean);
	//       if( CLHEP::RandFlat::shoot() > mean-multi ) multi++;
     }
#ifdef G4PHPDEBUG
     if( getenv("G4ParticleHPDebug") ) G4cout << "G4ParticleHPProduct::GetMultiplicity " << theMassCode << " " << theMass << " multi " << multi << " mean " << mean << G4endl;
#endif
  }

  fCache.Get().theCurrentMultiplicity = static_cast<G4int>(mean);

  return multi;
}


G4ReactionProductVector * G4ParticleHPProduct::Sample(G4double anEnergy, G4int multi)
{
  if(theDist == 0) { return 0; }
  G4ReactionProductVector * result = new G4ReactionProductVector;

  theDist->SetTarget(fCache.Get().theTarget);
  theDist->SetProjectileRP(fCache.Get().theProjectileRP);
  G4int i;
//  G4double eMax = GetTarget()->GetMass()+GetNeutron()->GetMass()
//                  - theActualStateQValue;
  G4ReactionProduct * tmp;
  theDist->ClearHistories();

  for(i=0;i<multi;i++)
  {
#ifdef G4PHPDEBUG
 if( getenv("G4PHPTEST") )
    if( getenv("G4ParticleHPDebug") && tmp != 0 )    G4cout << multi << " " << i << " @@@ G4ParticleHPProduct::Sample " << anEnergy << " Mass " << theMassCode << " " << theMass << G4endl;
#endif
    tmp = theDist->Sample(anEnergy, theMassCode, theMass);
    if(tmp != 0) { result->push_back(tmp); }
#ifndef G4PHPDEBUG //GDEB
    if( getenv("G4ParticleHPDebug") && tmp != 0 )   G4cout << multi << " " << i << " @@@ G4ParticleHPProduct::Sample " << tmp->GetDefinition()->GetParticleName() << " E= " << tmp->GetKineticEnergy() << G4endl; 
#endif
  }
  if(multi == 0) 
    {
      tmp = theDist->Sample(anEnergy, theMassCode, theMass);
      delete  tmp;
    }
  /*
  //080901 TK Comment out, too many secondaries are produced in deuteron reactions
  if(theTarget->GetMass()<2*GeV) // @@@ take care of residuals in all cases
  {
  tmp = theDist->Sample(anEnergy, theMassCode, theMass);
  tmp->SetDefinition(G4Proton::Proton());
  if(tmp != 0) { result->push_back(tmp); }
  }
  */

  return result;
}
