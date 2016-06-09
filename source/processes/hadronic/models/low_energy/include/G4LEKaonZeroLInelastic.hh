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
// $Id: G4LEKaonZeroLInelastic.hh,v 1.9 2003/07/01 15:49:03 hpw Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
//
// G4 Gheisha High Energy model class -- header file
// H. Fesefeldt, RWTH Aachen 23-October-1996
// Last modified: 10-December-1996

// A prototype of the Gheisha High Energy collision model.

#ifndef G4LEKaonZeroLInelastic_h
#define G4LEKaonZeroLInelastic_h 1

#include "G4LEKaonZeroInelastic.hh"
#include "G4LEAntiKaonZeroInelastic.hh"
#include "Randomize.hh"

class G4LEKaonZeroLInelastic : public G4InelasticInteraction  
{
  public: 
    G4LEKaonZeroLInelastic() 
    {
      SetMinEnergy( 0.0 );
      SetMaxEnergy( 25.*GeV );
    }

    virtual ~G4LEKaonZeroLInelastic(){ }

    G4HadFinalState * ApplyYourself(const G4HadProjectile &aTrack, G4Nucleus &targetNucleus )
    {
      if(G4UniformRand() < 0.50)
      {
         return theKaonZeroInelastic.ApplyYourself(aTrack, targetNucleus);
      }
      else
      {
         return theAntiKaonZeroInelastic.ApplyYourself(aTrack, targetNucleus);
      }
    } 
        
    G4LEKaonZeroInelastic theKaonZeroInelastic;
    G4LEAntiKaonZeroInelastic theAntiKaonZeroInelastic;
	

};
#endif                     
                                         

