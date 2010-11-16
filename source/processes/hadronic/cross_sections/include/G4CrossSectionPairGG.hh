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
// $Id: G4CrossSectionPairGG.hh,v 1.1 2010-11-16 13:38:15 gunter Exp $
// $ GEANT4 tag $Name: not supported by cvs2svn $
//
//   Class G4CrossSectionPairGG
//
//     Extend a cross section to higher energies using
//       G4GlauberGribovCrossSection at high energies.
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
#include "G4GlauberGribovCrossSection.hh"
#include <valarray>

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

  G4bool IsApplicable(const G4DynamicParticle* particle, const G4Element* element)
  { return IsZAApplicable(particle, element->GetZ(), element->GetN()); }  

  G4bool IsZAApplicable(const G4DynamicParticle* particle,
                        G4double ZZ, G4double AA);

  G4double GetCrossSection(const G4DynamicParticle* particle, 
                           const G4Element* element,
                           G4double temperature)
         { return GetIsoZACrossSection(particle, element->GetZ(), 
                                element->GetN(), temperature);    }

  G4double GetIsoZACrossSection(const G4DynamicParticle* particle,
                                G4double ZZ, G4double AA,
                                G4double /*aTemperature*/);



  void BuildPhysicsTable(const G4ParticleDefinition&);
  void DumpPhysicsTable(const G4ParticleDefinition&);
  
  private:
    G4VCrossSectionDataSet * theLowX;   
//    G4VCrossSectionDataSet * theHighX;
    G4GlauberGribovCrossSection * theHighX;
    G4double ETransition;
    typedef std::valarray<G4double> XS_factors;
    typedef std::pair<G4ParticleDefinition *, XS_factors > ParticleXScale;
    std::vector<ParticleXScale> scale_factors;

};

#endif
