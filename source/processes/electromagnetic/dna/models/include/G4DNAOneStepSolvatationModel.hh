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
// $Id: G4DNAOneStepSolvatationModel.hh 85244 2014-10-27 08:24:13Z gcosmo $
//
// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 


#ifndef G4DNAOneStepSolvatationModel_
#define G4DNAOneStepSolvatationModel_

#include "G4VEmModel.hh"

/**
  * When an electron reaches the highest energy domain of G4DNAOneStepSolvatationModel,
  * it is then automatically converted into a solvated electron and displace from its original
  * position using a published thermalization statistic.
  *
  * Article: Jintana Meesungnoen, Jean-Paul Jay-Gerin, Abdelali Filali-Mouhim, and Samlee Mankhetkorn (2002)
  * Low-Energy Electron Penetration Range in Liquid Water.
  * Radiation Research: November 2002, Vol. 158, No. 5, pp. 657-660.
  */

class G4DNAOneStepSolvatationModel : public G4VEmModel
{
public :
    G4DNAOneStepSolvatationModel(const G4ParticleDefinition* p = 0,
                             const G4String& nam = "DNASancheSolvatationModel");
    virtual ~G4DNAOneStepSolvatationModel();

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
    G4DNAOneStepSolvatationModel & operator=(const  G4DNAOneStepSolvatationModel &right);
    G4DNAOneStepSolvatationModel(const  G4DNAOneStepSolvatationModel&);
};

inline void G4DNAOneStepSolvatationModel::SetVerbose(int flag)
{
    fVerboseLevel = flag;
}

#endif
