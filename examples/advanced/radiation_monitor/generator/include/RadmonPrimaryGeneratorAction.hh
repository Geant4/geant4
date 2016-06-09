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
// File name:     RadmonPrimaryGeneratorAction.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPrimaryGeneratorAction.hh,v 1.1.2.2 2006/06/29 16:15:48 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//
// Description:   Primary particles generator user action
//

#ifndef   RADMONPRIMARYGENERATORACTION_HH
 #define  RADMONPRIMARYGENERATORACTION_HH
 
 // Include files
 #include "RadmonVLayoutObserver.hh"
 #include "G4VUserPrimaryGeneratorAction.hh"
 #include "G4ParticleGun.hh"
 #include "G4String.hh"
 #include "globals.hh"
 #include <map>
 
 // Forward declarations
 class RadmonVGeneratorsFactory;
 class RadmonVGeneratorLayout;
 class RadmonVGenerator;

 class RadmonPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction, public RadmonVLayoutObserver
 {
  public:
                                                RadmonPrimaryGeneratorAction(RadmonVGeneratorLayout * layout, RadmonVGeneratorsFactory * factory);
                                               ~RadmonPrimaryGeneratorAction();

   virtual void                                 GeneratePrimaries(G4Event * anEvent);

   virtual void                                 OnLayoutChange();

  private:
   inline void                                  Update(void);
   inline void                                  CleanUp(void);
   inline static G4String                       MakeKey(const G4String & source, const G4String & algorithm);

  // Hidden constructors and operators
                                                RadmonPrimaryGeneratorAction();
                                                RadmonPrimaryGeneratorAction(const RadmonPrimaryGeneratorAction & copy);
   RadmonPrimaryGeneratorAction &               operator=(const RadmonPrimaryGeneratorAction & copy);

  // Private data types
   typedef std::map<G4String, RadmonVGenerator *> GeneratorsMap;

  // Private attributes
   G4ParticleGun                                particlesGun;
   RadmonVGeneratorLayout *                     generatorLayout;
   RadmonVGeneratorsFactory *                   generatorsFactory;
   G4bool                                       needUpdate;
   GeneratorsMap                                generatorsMap;
 };
#endif /* RADMONPRIMARYGENERATORACTION_HH */
