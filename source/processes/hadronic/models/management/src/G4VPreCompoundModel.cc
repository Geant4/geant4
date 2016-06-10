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
// $Id: G4VPreCompoundModel.cc 80131 2014-04-02 14:35:47Z gcosmo $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class file
//
//      History: first implementation 1998
//
//      V.Ivanchenko 03.01.2012
//          Added G4ExcitationHandler pointer to the constructor and cleanup
// -----------------------------------------------------------------------------


#include "G4VPreCompoundModel.hh"

G4VPreCompoundModel::G4VPreCompoundModel(G4ExcitationHandler* ptr,
                                         const G4String& modelName):
  G4HadronicInteraction(modelName), theExcitationHandler(ptr)
{}

G4VPreCompoundModel::~G4VPreCompoundModel()
{}
 
void G4VPreCompoundModel::DeExciteModelDescription(std::ostream& outFile) const
{
   outFile << "description of DeExcite() for model derived from G4VPrecompoundModel missing"
           << "\n";

}
