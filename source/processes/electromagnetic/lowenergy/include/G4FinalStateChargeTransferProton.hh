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
// $Id: G4FinalStateChargeTransferProton.hh,v 1.2 2008-03-25 16:00:20 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// Final state policy class for proton charge transfer
// Reference for policy-based design:
// S. Chauvie et al., Geant4 physics processes for microdosimetry simulation:
// design foundation and implementation of the first set of models,
// IEEE Trans. Nucl. Sci., vol. 54, no. 6, Dec. 2007.
// Reference for implementation:
// http://www-pub.iaea.org/MTCD/publications/PDF/APID-VOL10.pdf 
// Further documentation available from http://www.ge.infn.it/geant4/

// -------------------------------------------------------------------


#ifndef G4FINALSTATECHARGETRANSFERCH_HH
#define G4FINALSTATECHARGETRANSFERCH_HH 1
 
#include "globals.hh"
#include "G4FinalStateProduct.hh"
#include <map>

class G4Track;
class G4Step;

 class G4FinalStateChargeTransferProton
 {
 public:
   
   G4FinalStateChargeTransferProton();
   
   ~G4FinalStateChargeTransferProton();
   
   const G4FinalStateProduct& GenerateFinalState(const G4Track& track, const G4Step& step);
   
 private:
   
   // Copy constructor and assignment operator to be added here
   
   G4String name;  
   G4double lowEnergyLimit;
   G4double highEnergyLimit;
   G4FinalStateProduct product;
  
   G4double protonBinding;

   std::map<G4String,G4double,std::less<G4String> > ioniPotentialMap;
};


#endif
