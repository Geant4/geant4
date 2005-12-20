//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4DNAGenericIonsManager.hh,v 1.2 2005-12-20 13:46:32 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef   G4DNAGENERICIONSMANAGER_HH
 #define  G4DNAGENERICIONSMANAGER_HH 1
 
 #include "globals.hh"
 
 #include <map>
 
 class G4Ions;
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
