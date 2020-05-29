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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4LENDGammaModel.hh                                               //
//  Date:   30 March 2020                                                     //
//  Author: Dennis H. Wright                                                  //
//                                                                            //
//  Description: model for inelastic scattering of gammas from nuclei         //
//               including gamma-induced fission.  This model is very similar //
//               to G4LENDCombinedModel except that it does not sample        //
//               elastic or capture reactions since there are no such data    //
//               for gammas in GND.                                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef G4LENDGammaModel_h
#define G4LENDGammaModel_h 1

#include "G4LENDModel.hh"

class G4LENDGammaCrossSection;

class G4LENDInelastic;
class G4LENDFission;

class G4LENDGammaModel : public G4LENDModel
{
  public: 
    G4LENDGammaModel(G4ParticleDefinition* pd);
    ~G4LENDGammaModel(){;};

    void BuildPhysicsTable(const G4ParticleDefinition&);

    G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
                                   G4Nucleus& aTargetNucleus);

    G4bool HasData(const G4DynamicParticle*, G4int iZ, G4int iA, G4int iM, 
                   const G4Isotope*, const G4Element*, const G4Material*);
     
  private: 
    G4LENDGammaCrossSection* crossSection;
    G4LENDInelastic* inelastic;
    G4LENDFission* fission;
    G4LENDModel* channels[2];
};

#endif
