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
// Class G4TCatalog
//
// Class description:
//
// Handles the management of the catalog.txt file (loads the items or
// inserts new ones).
//
// History:
// Created by Roman Atachiants, 18/08/2009
// Modified:
// M.Kosov 25/05/2010: catalog reading/wrighting correction
//
// --------------------------------------------------------------------
#ifndef G4TCatalog_H_
#define G4TCatalog_H_

#include "../CommonHeaders.h"

class G4TCatalog : public TObject
{
  private:

  TString fFileName;

  public:

  G4TCatalog() : fFileName("./database/catalog.txt")  { }
  virtual ~G4TCatalog () {}

  vector<TString>   Load();
  Bool_t    ContainsLine(TString const& line);
  void      Insert(TString const& name);

  vector<TString>   GetPublications();

  ClassDef(G4TCatalog, 1)  //The class for the data catalog
};

R__EXTERN G4TCatalog * gCatalog;
 
#endif
