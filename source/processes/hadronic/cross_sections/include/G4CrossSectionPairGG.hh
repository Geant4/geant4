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
// $Id: G4CrossSectionPairGG.hh 93682 2015-10-28 10:09:49Z gcosmo $
// $ GEANT4 tag $Name: not supported by cvs2svn $
//
//   Class G4CrossSectionPairGG
//
//     Extend a cross section to higher energies using
//       G4ComponentGGHadronNucleusXsc at high energies.
//       Smoothly join cross section sets by scaling GG at a given 
//       transition energy to match the given low energy cross section.
//
//  Author:  Gunter Folger
//           November 2010
//
#ifndef G4CrossSectionPairGG_h
#define G4CrossSectionPairGG_h

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ParticleDefinition.hh"
#include <valarray>
#include <iostream>

class G4NistManager;
class G4ComponentGGHadronNucleusXsc;

class G4CrossSectionPairGG : public G4VCrossSectionDataSet
{
  private:
   G4CrossSectionPairGG();
   G4CrossSectionPairGG(const G4CrossSectionPairGG&);
   G4CrossSectionPairGG& operator=(const G4CrossSectionPairGG&);

  public:
  
  G4CrossSectionPairGG(G4VCrossSectionDataSet * low,
//  	             G4VCrossSectionDataSet * high,
		     G4double Etransit);
			      
  virtual ~G4CrossSectionPairGG();

  virtual void CrossSectionDescription(std::ostream&) const;

  virtual
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z, 
			     const G4Material* mat = 0);

  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, G4int Z,
				  const G4Material* mat = 0);

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);
  virtual void DumpPhysicsTable(const G4ParticleDefinition&);
  
  private:

    G4NistManager* NistMan;
    G4VCrossSectionDataSet * theLowX;   
    G4ComponentGGHadronNucleusXsc * theHighX;
    G4double ETransition;
    typedef std::valarray<G4double> XS_factors;
    typedef std::pair<const G4ParticleDefinition *, XS_factors > ParticleXScale;
    std::vector<ParticleXScale> scale_factors;

};

#endif
