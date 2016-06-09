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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4LElasticB.hh,v 1.1 2005/12/14 18:13:06 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// G4 Model: Low energy elastic scattering with 4-momentum balance 
// Derived fron G4LElastic of F.W. Jones, TRIUMF, 04-JUN-96
//  
// Modified:
// 14-Dec-05 V.Ivanchenko rename the class
//
// use -scheme for elastic scattering: HPW, 20th June 1997
// most of the code comes from the old Low-energy Elastic class
//


#ifndef G4LElasticB_h
#define G4LElasticB_h 1
 
// Class Description
// Final state production model for hadron nuclear elastic scattering; 
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

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


class G4LElasticB : public G4HadronicInteraction
{
public:

   G4LElasticB() : G4HadronicInteraction()
   {
      SetMinEnergy( 0.0*GeV );
      SetMaxEnergy( DBL_MAX );
   }

   ~G4LElasticB() {};
 
   G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, 
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
