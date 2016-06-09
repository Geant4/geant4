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
// $Id: G4ParallelGeometrySampler.hh,v 1.8 2006/06/29 21:10:10 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// Class G4ParallelGeometrySampler
//
// Class description:
//
// Class description:
// This class inherits from G4VSampler. It is used for scoring and 
// importance sampling in a parallel geometry.
// See also the description in G4VSampler.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelGeometrySampler_hh
#define G4ParallelGeometrySampler_hh G4ParallelGeometrySampler_hh

#include "G4Types.hh"
#include "G4VSampler.hh"
#include "G4ParallelWorld.hh"
#include "G4VSamplerConfigurator.hh"

class G4ParallelTransportConfigurator;
class G4PScoreConfigurator;
class G4PImportanceConfigurator;
class G4WeightCutOffConfigurator;
class G4VGCellFinder;
class G4PWeightWindowConfigurator;
class G4WeightWindowStore;

class G4ParallelGeometrySampler : public G4VSampler
{

public:  // with description
 
  G4ParallelGeometrySampler(G4VPhysicalVolume &worldvolume,
                            const G4String &particlename);
  virtual ~G4ParallelGeometrySampler();

  virtual void PrepareScoring(G4VScorer *Scorer);
  virtual void PrepareImportanceSampling(G4VIStore *istore,
                                 const G4VImportanceAlgorithm 
                                 *ialg);
  virtual void PrepareWeightRoulett(G4double wsurvive, 
                            G4double wlimit,
                            G4double isource);
  
  virtual void PrepareWeightWindow(G4VWeightWindowStore *wwstore,
                                   G4VWeightWindowAlgorithm *wwAlg,
                                   G4PlaceOfAction placeOfAction);


  virtual void Configure();

  virtual void ClearSampling();
  virtual G4bool IsConfigured() const;
  
private:

  G4ParallelGeometrySampler(const G4ParallelGeometrySampler &);
  G4ParallelGeometrySampler &
  operator=(const G4ParallelGeometrySampler &);

private:

  G4String fParticleName;
  G4ParallelWorld fParallelWorld;
  G4ParallelTransportConfigurator *fParallelTransportConfigurator;
  G4PImportanceConfigurator *fPImportanceConfigurator;
  G4PScoreConfigurator *fPScoreConfigurator;
  G4VGCellFinder *fGCellFinder;
  G4WeightCutOffConfigurator *fWeightCutOffConfigurator;
  G4VIStore *fIStore;
  G4PWeightWindowConfigurator *fPWeightWindowConfigurator;
  G4VWeightWindowStore *fWWStore;
  G4bool fIsConfigured;
  G4Configurators fConfigurators;
};
  
#endif
