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
// $Id: G4CrossSectionChargeTransferExp.hh,v 1.2 2008-03-25 16:00:20 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// Date         Name              Modification
// 28 Dec 2007  M.G. Pia          Created
//
// -------------------------------------------------------------------

// Class description:

// Total cross section for charge transfer from experimental data interpolation
//            http://www-cfadc.phy.ornl.gov/astro/ps/data/
// Further documentation available from http://www.ge.infn.it/geant4/

// -------------------------------------------------------------------

#ifndef G4CROSSSECTIONCHARGETRANSFEREXP_HH
#define G4CROSSSECTIONCHARGETRANSFEREXP_HH 1
 
#include "globals.hh"
#include <map>

class G4Track;
class G4VEMDataSet;
 
class G4CrossSectionChargeTransferExp
{
public:
  
  G4CrossSectionChargeTransferExp();
  
  virtual ~G4CrossSectionChargeTransferExp();
  
  G4double CrossSection(const G4Track&);
  
  // Copy constructor and assignment operator to be added here
    
private:
   
  G4String name;  
  G4double lowEnergyLimit;
  G4double highEnergyLimit;
  
  G4double unit1;
  G4double unit2;
 
  std::map<G4String,G4VEMDataSet*,std::less<G4String> > crossMap;

  G4VEMDataSet* LoadData(const G4String& dataFile);

};

#endif
