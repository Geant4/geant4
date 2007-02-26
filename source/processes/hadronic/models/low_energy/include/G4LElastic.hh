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
// $Id: G4LElastic.hh,v 1.12 2007-02-26 19:05:01 dennis Exp $
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
#include "G4LightMedia.hh"
#include "G4Step.hh"
#include "G4TrackStatus.hh"
#include "G4HadronicInteraction.hh"


class G4LElastic : public G4HadronicInteraction
{
public:

   G4LElastic() : G4HadronicInteraction("G4LElastic")
   {
      SetMinEnergy( 0.0*GeV );
      SetMaxEnergy( DBL_MAX );
   }

   ~G4LElastic() {};
 
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
