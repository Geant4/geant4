// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LElastic.hh,v 1.3 1999-11-16 16:30:05 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Model: Low energy elastic scattering -- header file
// F.W. Jones, TRIUMF, 04-JUN-96
//  
// For further comments see G4LElastic.cc.
//
// use -scheme for elastic scattering: HPW, 20th June 1997
// most of the code comes from the old Low-energy Elastic class
//


#ifndef G4LElastic_h
#define G4LElastic_h 1
 
#include "globals.hh"
#include "Randomize.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4LPhysicsFreeVector.hh"
//#include "G4HadronicCrossSections.hh"
#include "G4LightMedia.hh"
#include "G4Step.hh"
#include "G4TrackStatus.hh"
#include "G4HadronicInteraction.hh"


class G4LElastic : public G4HadronicInteraction
{
public:

   G4LElastic() : G4HadronicInteraction()
   {
      SetMinEnergy( 0.0*GeV );
      SetMaxEnergy( DBL_MAX );
   }

   ~G4LElastic() {};
 
   G4VParticleChange* ApplyYourself(const G4Track& aTrack, 
                                    G4Nucleus & targetNucleus);

private:


   G4LightMedia LightMedia;

   G4int Rtmi(G4double* x, G4double xli, G4double xri, G4double eps, 
              G4int iend,
              G4double aa, G4double bb, G4double cc, G4double dd, 
              G4double rr);

   G4double Fctcos(G4double t, 
                   G4double aa, G4double bb, G4double cc, G4double dd, 
                   G4double rr);

   void Defs1(G4double p, G4double px, G4double py, G4double pz, 
              G4double pxinc, G4double pyinc, G4double pzinc, 
              G4double* pxnew, G4double* pynew, G4double* pznew);
};
#endif
