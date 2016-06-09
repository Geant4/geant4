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
// $Id: G4DNASancheSolvatationModel.hh 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4SancheSolvatationModel_
#define G4SancheSolvatationModel_

#include "G4VEmModel.hh"

/**
  * When an electron reaches the highest energy domain of G4DNASancheSolvatationModel,
  * it is then automatically converted into a solvated electron and displace from its original
  * position using a published thermalization statistic.
  *
  */

class G4DNASancheSolvatationModel : public G4VEmModel
{
public :
    G4DNASancheSolvatationModel(const G4ParticleDefinition* p = 0,
                             const G4String& nam = "DNASancheSolvatationModel");
    virtual ~G4DNASancheSolvatationModel();

    virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

    virtual G4double CrossSectionPerVolume(  const G4Material* material,
            const G4ParticleDefinition* p,
            G4double ekin,
            G4double emin,
            G4double emax);

    virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                   const G4MaterialCutsCouple*,
                                   const G4DynamicParticle*,
                                   G4double tmin,
                                   G4double maxEnergy);

    inline void SetVerbose(int);

protected:
    // Water density table
    const std::vector<G4double>* fpWaterDensity;

    G4ThreeVector RadialDistributionOfProducts(G4double Rrms) const ;
    G4ParticleChangeForGamma* fParticleChangeForGamma;

    G4bool fIsInitialised;
    G4int fVerboseLevel;

private :
    G4DNASancheSolvatationModel & operator=(const  G4DNASancheSolvatationModel &right);
    G4DNASancheSolvatationModel(const  G4DNASancheSolvatationModel&);
};

inline void G4DNASancheSolvatationModel::SetVerbose(int flag)
{
    fVerboseLevel = flag;
}

#endif
