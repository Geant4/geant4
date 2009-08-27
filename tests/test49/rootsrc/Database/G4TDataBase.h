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
// Class G4TDataBase
//
// Class description:
//
// The class to handle loading and saving of the G4TData objects, it
// is recommended to only use this class and not use directly the Save()
// and Load() methods of G4TData class.
//
// History:
// Roman Atachiants, 18/08/2009 - initial version
//
// --------------------------------------------------------------------

#ifndef G4TDataBase_H_
#define G4TDataBase_H_

#include "../CommonHeaders.h"
#include "G4TCatalog.h"
#include "G4TData.h"


class G4TDataBase : public TObject {

  private:
	  TString fDirectory;

  public:

	  G4TDataBase() : fDirectory("./database/") {  }
	  virtual ~G4TDataBase () {}

	  // Access Methods
	  G4TData*		LoadData(G4TData* object, Int_t secondaryPDG = 0);
	  G4TData*		LoadData(TString const& filename, Int_t secondaryPDG = 0);
	  void			SaveData(G4TData* object);

	  // Getters/Setters
	  TString		GetDirectory() const;
	  void			SetDirectory(TString fDirectory);


	  ClassDef(G4TDataBase, 1)  //The class for Geant4 Testing Database DAL
};


R__EXTERN G4TDataBase *gTestDB;

#endif




