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
// $Id: $
//
// Author:  D.H. Wright
// Date:    2 February 2011
//
// Description: model of muon nuclear interaction in which a gamma from
//              the virtual photon spectrum interacts in the nucleus as
//              a real gamma at low energies and as a pi0 at high energies.
//              Kokoulin's muon cross section and equivalent gamma spectrum
//              are used.
//
 
#ifndef G4MuonVDNuclearModel_h
#define G4MuonVDNuclearModel_h 1 

#include "G4HadronicInteraction.hh"

class G4CascadeInterface;
class G4TheoFSGenerator; 
class G4LundStringFragmentation;
class G4ExcitedStringDecay;
class G4KokoulinMuonNuclearXS;
class G4ElementData; 

class G4MuonVDNuclearModel : public G4HadronicInteraction
{
public:
    
  explicit G4MuonVDNuclearModel();
  virtual ~G4MuonVDNuclearModel();

  virtual G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
					 G4Nucleus& targetNucleus);

  virtual void ModelDescription(std::ostream& outFile) const;

private:
    
  G4DynamicParticle* CalculateEMVertex(const G4HadProjectile& aTrack,
				       G4Nucleus& targetNucleus);

  void CalculateHadronicVertex(G4DynamicParticle* incident,
			       G4Nucleus& target);

  void MakeSamplingTable();

  G4MuonVDNuclearModel & operator=(const G4MuonVDNuclearModel &right) = delete;
  G4MuonVDNuclearModel(const  G4MuonVDNuclearModel&) = delete;

  static const G4int NBIN = 800;
  static const G4int nzdat = 5;
  static const G4int ntdat = 73;

  static const G4int zdat[nzdat];
  static const G4double adat[nzdat];
  static const G4double tdat[ntdat];
  
  static G4ElementData* fElementData;             

  G4double CutFixed;
  G4bool isMaster;

  G4KokoulinMuonNuclearXS* muNucXS;

  G4TheoFSGenerator* ftfp;
  G4LundStringFragmentation* theFragmentation;
  G4ExcitedStringDecay* theStringDecay;
  G4CascadeInterface* bert;
};
 
#endif

