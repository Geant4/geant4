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
// -------------------------------------------------------------------
//      GEANT 4 test file  
// 
//
//      File name:   ParticleScattererTest.cc   
//
//      Authors:       Yuri Smirnov (Youri.Smirnov@cern.ch
//                                   smirnov@lcta6.jinr.dubna.su
//                                    smirnovy@cv.jinr.dubna.su)
//
//      Creation date:       June 1998
//
//      Last Modification:  13 July 1998 
//   
//      Version:            beta       
// ------------------------------------------------------------------- 

#include <iostream>    
#include <stdlib.h>     
#include <math.h> 
#include "G4ParticleScatterer.hh" 
#include "G4CrossSectionBase.hh"
#include "globals.hh"
#include "G4ElasticScatterer.hh"
#include "G4InElasticScatterer.hh"
#include "G4ElasticScatterer.hh"
#include "G4Annihilator.hh"

#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh" 

int main()
{
  
 G4KineticTrack Projectile;  
 G4KineticTrack Target; 
 G4KineticTrackVector* theOutputTrackVector;
  
//CMS test: working 
// G4LorentzVector  Mo1 (0.0*GeV, 0.0*GeV,  0.1*GeV, 0.943586*GeV);
// G4LorentzVector  Mo2 (0.0*GeV, 0.0*GeV, -0.1*GeV, 0.944872*GeV);
// G4ThreeVector    Po1 (0.0*fermi,  0.0*fermi,   0.1*fermi);
// G4ThreeVector    Po2 (0.0*fermi,  0.0*fermi,  -0.1*fermi);

// !!!  Working Tested version !!!   
//G4LorentzVector  Mo1 (0.5*GeV, 0.0*GeV,  0.1*GeV, 1.1403543*GeV);
//G4LorentzVector  Mo2 (0.0*GeV, 0.0*GeV, -0.1*GeV, 0.9448725*GeV);
//G4ThreeVector    Po1 (0.0*fermi,  0.0*fermi,   0.1*fermi);
//G4ThreeVector    Po2 (0.0*fermi,  0.0*fermi,  -0.1*fermi);

G4LorentzVector  Mo1 (0.0*GeV, 0.0*GeV,  0.1*GeV, 0.943586*GeV);
G4LorentzVector  Mo2 (0.0*GeV, 0.0*GeV, -0.1*GeV, 0.9448725*GeV);
G4ThreeVector    Po1 (1.0*fermi,  1.0*fermi,   1.0*fermi);
G4ThreeVector    Po2 (0.5*fermi, -1.0*fermi,   3.0*fermi);
cout << "G4LorentzVector1 = " << Mo1 <<  G4endl;
cout << "G4LorentzVector2 = " << Mo2 <<  G4endl;    
 
 
G4ParticleDefinition* Defin1 = G4Proton::ProtonDefinition();
G4ParticleDefinition* Defin2 = G4Neutron::NeutronDefinition();   

 Projectile.SetDefinition(Defin1);
 Target.SetDefinition(Defin2);

 Projectile.Set4Momentum(Mo1); 
 Target.Set4Momentum(Mo2); 
 Projectile.SetPosition(Po1);
 Target.SetPosition(Po2);
 
 G4String Nam1 = Defin1->GetParticleName();
 G4String Nam2 = Defin2->GetParticleName();   

cout << "======================================================= " << G4endl;
cout << " ParticleScattererTest: " << G4endl;
cout << "Input information:" << G4endl;
cout << "Particle1 = " << Nam1 << "  Particle2 = " << Nam2 << G4endl;
cout << "G4LorentzVector1 = " << Mo1 <<  G4endl;
cout << "G4LorentzVector2 = " << Mo2 <<  G4endl;
cout << "G4ThreeVector1 = " << Po1 <<  G4endl;
cout << "G4ThreeVector2 = " << Po2 <<  G4endl; 

  G4LorentzVector Moment1(Projectile.Get4Momentum());
  G4LorentzVector Moment2(Projectile.Get4Momentum());

  G4CrossSectionBase CSB;   
  G4ElasticScatterer ElS;   
  G4InElasticScatterer InS;   
 
  G4ParticleScatterer Scat1(&CSB, &ElS, &InS, NULL);  
  G4double MyTime1 = Scat1.GetTimeToInteraction(Projectile,Target);
                                          
cout << "======================================================= " << G4endl;
cout << " ParticleScattererTestPart:  After GetTimeToInteraction :" << G4endl;
cout << "Input information:" << G4endl;
cout << "Particle1 = " << Nam1 << " Particle2 = " << Nam2 << G4endl;
cout << "G4LorentzVector1 = " << Mo1 <<  G4endl;
cout << "G4LorentzVector2 = " << Mo2 <<  G4endl;
cout << "G4ThreeVector1 = " << Po1 <<  G4endl;
cout << "G4ThreeVector2 = " << Po2 <<  G4endl;
cout << "=== TimeToInteraction = " << MyTime1 << G4endl;

// G4bool IsScat = Scat1.IsScattering(Projectile, Target);

//cout << "======================================================= " << G4endl;
//cout << " ParticleScateerTest: After IsScattering :" << G4endl;
//cout << "Output information:" << G4endl;
//cout << "Particle1 = " << Nam1 << "Particle2 = " << Nam2 << G4endl;
//cout << "G4LorentzVector1 = " << Mo1 <<  G4endl;
//cout << "G4LorentzVector2 = " << Mo2 <<  G4endl;
//cout << "G4ThreeVector1 = " << Po1 <<  G4endl;
//cout << "G4ThreeVector2 = " << Po2 <<  G4endl; 
//cout << "======================================================= " << G4endl;
//cout << " IsScattering = " <<  IsScat << G4endl;
     
 theOutputTrackVector = Scat1.Scatter(Projectile, Target);
 
cout << "************************************* " << G4endl;
cout << " ParticleScattererTest: After Scatter  " << G4endl;
cout << " See theOutputTrackVector " << G4endl;
cout << "************************************* " << G4endl;
       G4double px, py, pz, e;
       px = py =pz = e  = 0;
        for(int c1 = 0; c1 < theOutputTrackVector->length(); c1++) 
              {  
                G4KineticTrack* Trk = theOutputTrackVector->at(c1);
                                Moment1=Trk->Get4Momentum(); 
                                Po1=Trk->GetPosition(); 
                                Defin1=Trk->GetDefinition();
                                Nam1 = Defin1->GetParticleName();
cout << " Current Particle Number =" << c1  << G4endl;
cout << " Current Particle Name : " << Nam1 << G4endl;
cout << " G4LorentzVector = " << Moment1 <<  G4endl;
cout << " G4ThreeVector = " << Po1 <<  G4endl;
              }   

  return 0;              
}

	
