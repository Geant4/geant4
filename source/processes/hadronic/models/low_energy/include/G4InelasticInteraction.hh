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
// $Id$
//
// Hadronic Process: Inelastic Interaction 
// This class is an abstract base class, since the pure virtual
// function ApplyYourself has not been defined yet.
// original by H.P. Wellisch
// Modified by J.L. Chuma, TRIUMF, 22-Nov-1996
// Modified by J.L. Chuma  27-Mar-1997
// Modified by J.L. Chuma  30-Apr-1997
// Modified by J.L. Chuma  05-Aug-1997  to pass the original incident particle to
//                                      CalculateMomenta
// Modified by J.L. Chuma  05-Jun-1998  to include quasiElastic flag to allow for
//                                      TwoBody to be called directly, bypassing
//                                      TwoCluster, and allowing TwoCluster to be
//                                      called with no secondaries
// 23-Jan-2009 V.Ivanchenko move constructor and destructor to the body
// 29-Aug-2009 V.Ivanchenko moveded G4ReactionDynamics instance from the based class 
 
#ifndef G4InelasticInteraction_h
#define G4InelasticInteraction_h 1

#include "globals.hh"
#include "G4FastVector.hh"
#include "G4HadronicInteraction.hh"
#include "G4ReactionProduct.hh"
#include "G4ParticleTypes.hh" 
#include "Randomize.hh"
#include "G4ReactionDynamics.hh"
#include "G4VIsotopeProduction.hh"

class G4IsoResult;
class G4IsoParticleChange;


class G4InelasticInteraction : public G4HadronicInteraction
{
  public:
    
    G4InelasticInteraction(const G4String& name = "LEInelastic");
    virtual ~G4InelasticInteraction();

    void RegisterIsotopeProductionModel(G4VIsotopeProduction* aModel)
      {theProductionModels.push_back(aModel);}

    void TurnOnIsotopeProduction() {isotopeProduction = true;}

    static G4IsoParticleChange* GetIsotopeProductionInfo(); 

    virtual const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const;

  protected:
    
    G4double Pmltpc(G4int np, G4int nm, G4int nz, G4int n,
                    G4double b, G4double c );
    
    G4bool MarkLeadingStrangeParticle(const G4ReactionProduct &currentParticle,
                                      const G4ReactionProduct &targetParticle,
                                      G4ReactionProduct &leadParticle);
    
    void SetUpPions(const G4int np, const G4int nm, const G4int nz,
                    G4FastVector<G4ReactionProduct,GHADLISTSIZE>& vec,
                    G4int& vecLen);
    
    void Rotate(G4FastVector<G4ReactionProduct,GHADLISTSIZE>& vec, G4int& vecLen);

    void GetNormalizationConstant(const G4double availableEnergy,
                                  G4double& n, G4double& anpn);
    
    void CalculateMomenta(G4FastVector<G4ReactionProduct,GHADLISTSIZE>& vec,
                          G4int& vecLen,
                          const G4HadProjectile* originalIncident,
                          const G4DynamicParticle* originalTarget,
                          G4ReactionProduct& modifiedOriginal,
                          G4Nucleus& targetNucleus,
                          G4ReactionProduct& currentParticle,
                          G4ReactionProduct& targetParticle,
                          G4bool& incidentHasChanged,
                          G4bool& targetHasChanged,
                          G4bool quasiElastic);
    
    void SetUpChange(G4FastVector<G4ReactionProduct,GHADLISTSIZE>& vec,
                     G4int& vecLen,
                     G4ReactionProduct& currentParticle,
                     G4ReactionProduct& targetParticle,
                     G4bool& incidentHasChanged);

    void DoIsotopeCounting(const G4HadProjectile* theProjectile,
                           const G4Nucleus& aNucleus);

    G4IsoResult* ExtractResidualNucleus(const G4Nucleus& aNucleus);
    G4bool isotopeProduction;

    G4ReactionDynamics theReactionDynamics;

  private:

    G4double cache;
    G4ThreeVector what;

    std::vector<G4VIsotopeProduction*> theProductionModels;
    static G4IsoParticleChange* theIsoResult;
    static G4IsoParticleChange* theOldIsoResult;
};
 
#endif
 
