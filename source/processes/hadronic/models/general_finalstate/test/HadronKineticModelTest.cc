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
//
// $Id: HadronKineticModelTest.cc,v 1.1 2003-10-08 12:04:17 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ios.hh"

#include "g4templates.hh"

#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

#include "G4HadronKineticModel.hh"
#include "G4Proton.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4Fancy3DNucleus.hh"



int main()
{
   G4cout << G4endl << "The new test program" << G4endl;

   G4HadronKineticModel anHadronKineticModel;



//    
//    Define the incident particle and the target nucleus (p + C)
//    

   G4ParticleDefinition* theProtonDefinition = G4Proton::ProtonDefinition();
   G4double              theProtonFormationTime = 0.0*ns;
   G4ThreeVector         theProtonInitialPosition 
                         (0.0*cm, 0.0*cm, 0.0*cm);
   G4LorentzVector       theProton4Momentum 
                         (0.0*MeV, 0.0*MeV, 1.0*MeV, 938.27*MeV);  
   G4KineticTrack        aProton(theProtonDefinition, 
                                 theProtonFormationTime,
                                 theProtonInitialPosition,
                                 theProton4Momentum);
                          
   G4Fancy3DNucleus      theCarbon;
   G4double              theA = 12.0;
   G4double              theZ =  6.0;
   theCarbon.Init(theA, theZ);
   
   G4cout << G4endl << " Nucleus mass number = " << theCarbon.GetMassNumber() << G4endl
                << " Nucleus mass = " << theCarbon.GetMass() << G4endl
                << " Nucleus charge = " << theCarbon.GetCharge() << G4endl;



                
//    
//    Setup the lits of projectile particles
//    

   G4KineticTrackVector theProjectileParticles;
   theProjectileParticles.insert(new G4KineticTrack(aProton));
//   theProjectileParticles.insert(new G4KineticTrack(aProton));
   
   
   
//    
//    Initialize some G4HadronKineticModel data members
//    

   anHadronKineticModel.SetProjectileTracks(theProjectileParticles);
   anHadronKineticModel.Set3DNucleus(&theCarbon);
   
   
   
//    
//    Do a time Step
//    

   anHadronKineticModel.DoTimeStep();
}    
    
