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
// $Id: G4DNAGenericIonsManager.hh 73124 2013-08-19 07:53:33Z gcosmo $

#ifndef   G4DNAGENERICIONSMANAGER_HH
 #define  G4DNAGENERICIONSMANAGER_HH 1
 
 #include "globals.hh"
 
 #include <map>
 
 class G4DNAIons;
 class G4ParticleDefinition;
 
 class G4DNAGenericIonsManager
 {
  public:
   static G4DNAGenericIonsManager *      Instance(void);
   G4ParticleDefinition *                GetIon(const G4String & name);

  private:
                                         G4DNAGenericIonsManager();
                                        ~G4DNAGenericIonsManager();

                                         G4DNAGenericIonsManager(const G4DNAGenericIonsManager &);
   const G4DNAGenericIonsManager        &operator=(const G4DNAGenericIonsManager &);

   static G4DNAGenericIonsManager *      theInstance;
   
   typedef std::map<G4String, G4ParticleDefinition *> IonsMap;

   IonsMap                               map;
 };
#endif /* G4DNAGENERICIONSMANAGER_HH */
