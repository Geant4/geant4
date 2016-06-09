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
// $Id: G4CrossSectionIonisationRuddPartial.hh,v 1.1 2007/11/08 21:35:31 pia Exp $
// GEANT4 tag $Name: geant4-09-01 $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// Date         Name              Modification
// 28 Apr 2007  M.G. Pia          Created in compliance with design described in TNS paper
//
// -------------------------------------------------------------------

// Class description:
// Geant4-DNA Cross total cross section for ionisation cross section in water
// Reference: TNS Geant4-DNA paper
// S. Chauvie et al., Geant4 physics processes for microdosimetry simulation:
// design foundation and implementation of the first set of models,
// IEEE Trans. Nucl. Sci., vol. 54, no. 6, Dec. 2007.
// Reference for implementation model: NIM. 155, pp. 145-156, 1978
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------


#ifndef G4CROSSSECTIONIONISATIONRuddPARTIAL_HH
#define G4CROSSSECTIONIONISATIONRuddPARTIAL_HH 1
 
#include "globals.hh"
#include <map>
#include <functional>
#include "G4DNACrossSectionDataSet.hh"
 
class G4Track;

class G4CrossSectionIonisationRuddPartial
{
public:
  
  G4CrossSectionIonisationRuddPartial();
  
  ~G4CrossSectionIonisationRuddPartial();
  
  // Partial cross section 
  // G4double CrossSection(G4double energy, const G4String& particle);
  G4double CrossSection(const G4Track& track);
			
  // Sum of partial cross sections at a given energy value for a particle type
  G4double Sum(G4double energy, const G4String& particle);

  G4int RandomSelect(G4double energy,const G4String& particle );
  
  // Copy constructor and assignment operator to be added here
    
private:
   
  G4String name;  
  G4double lowEnergyLimitDefault;
  G4double highEnergyLimitDefault;

  std::map<G4String,G4double,std::less<G4String> > lowEnergyLimit;
  std::map<G4String,G4double,std::less<G4String> > highEnergyLimit;

  typedef std::map<G4String,G4String,std::less<G4String> > MapFile;
  MapFile tableFile;

  typedef std::map<G4String,G4DNACrossSectionDataSet*,std::less<G4String> > MapData;
  MapData tableData;



  // G4DNACrossSectionDataSet* table;
 
};

#endif
