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
// $Id: G4RPGKLongInelastic.hh 79697 2014-03-12 13:10:09Z gcosmo $
//
// Author: D. H. Wright
// Date:   18 June 2007
//

// Class Description:
// Final state production model for K0L inelastic scattering
// using the re-parameterized Gheisha model.


#ifndef G4RPGKLongInelastic_h
#define G4RPGKLongInelastic_h 1

#include <CLHEP/Units/SystemOfUnits.h>
#include "G4RPGInelastic.hh"
#include "G4RPGKZeroInelastic.hh"
#include "G4RPGAntiKZeroInelastic.hh"
#include "Randomize.hh"


class G4RPGKLongInelastic : public G4RPGInelastic
{
  public: 
    G4RPGKLongInelastic() : G4RPGInelastic("G4RPGKLongInelastic")  
    {
      SetMinEnergy( 0.0 );
      SetMaxEnergy( 25.*CLHEP::GeV );
    }

    virtual ~G4RPGKLongInelastic(){ }

    G4HadFinalState * ApplyYourself(const G4HadProjectile &aTrack, G4Nucleus &targetNucleus )
    {
      if(G4UniformRand() < 0.50)
      {
         return theKZeroInelastic.ApplyYourself(aTrack, targetNucleus);
      }
      else
      {
         return theAntiKZeroInelastic.ApplyYourself(aTrack, targetNucleus);
      }
    } 
        
    G4RPGKZeroInelastic theKZeroInelastic;
    G4RPGAntiKZeroInelastic theAntiKZeroInelastic;
	

};
#endif                     
                                         

