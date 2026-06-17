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
// G4PhysicsModelCatalog
//
// Class description:
//
// Singleton, collection of physics models, to be used by models and G4Track.
//
// Author: M.Asai (SLAC), 26 September 2013
//
// Revised in August 2021 and May 2026 by A.Ribon (CERN).
// --------------------------------------------------------------------
#ifndef G4PHYSICSMODELCATALOG_HH
#define G4PHYSICSMODELCATALOG_HH

#include <vector>

#include "G4String.hh"
#include "globals.hh"


class G4PhysicsModelCatalog {
  public:
    static void Initialize();
    ~G4PhysicsModelCatalog()                              = default;
    G4PhysicsModelCatalog( const G4PhysicsModelCatalog& ) = delete;
    G4PhysicsModelCatalog& operator=( const G4PhysicsModelCatalog& ) = delete;

    static const G4String GetModelNameFromID( const G4int modelID );
    static const G4String GetModelNameFromIndex( const G4int modelIndex );
    static G4int GetModelID( const G4int modelIndex );
    static G4int GetModelID( const G4String& modelName );
    static G4int GetModelIndex( const G4int modelID );
    static G4int GetModelIndex( const G4String& modelName );
    // For Geant4 native (i.e. not custom/user-defined) physics models
    // there are two integer values: the model ID and the model Index.
    // The model ID is a unique integer which identifies not only the
    // Geant4 native physics model, but also the category to which it
    // belongs to; future Geant4 native physics models should be given
    // a model ID consistent with their category.
    // In the Geant4 code, it should always be the model ID, not the
    // model Index, which is used.
    // The model Index is the index of the vector of either model IDs, or
    // model names (these two vectors have the same size).
    // The model Index for a Geant4 native physics model does not have any
    // meaning in itself: it depends only on the order in which the vectors
    // are filled.
    // The model Index is useful for plotting because the index of the vector
    // has contiguous, small non-negative integer values, whereas the modelID
    // has non-contiguous, large, positive integer values (which is
    // unconvenient for plotting).
    // The idea is that, starting from Geant4 version 11.0, all the three
    // identifiers (modelID, index, name) for Geant4 native physics models
    // remain the same regardless of the physics list, application, and
    // version of Geant4.

    static G4int Entries();
    // The size of the two vectors (of model IDs and model names - for Geant4
    // native physics models) are required to be the same.

    static void PrintAllInformation();
    // Print all information of this class about all models - both Geant4
    // native and custom/user-defined physics models - i.e. each entry of
    // the two vectors - the one of model-IDs and the one of model-names -
    // for Geant4 native physics models, and the vector of model-names for
    // the custom/user-defined physics models.

    static G4int GetMinAllowedModelIDValue();
    static G4int GetMaxAllowedModelIDValue();
    // Returns the two limits, min and max respectively, that the modelID value
    // can have for Geant4 native physics models.

    static G4int RegisterCustomModel( const G4String& modelName );  
    // For custom/user-defined physics models, the user is responsible to call
    // this method during initialisation to register them, one by one, by
    // providing their names; the return value of the method is the unique
    // modelID assigned to the physics model whose name is specified in the
    // argument of the method.
    // Notes:
    // - Custom/user-defined physics models have modelID larger than
    //   GetMaxAllowedModelIDValue();
    // - Only the names of custom/user-defined physics models are kept in a
    //   separate vector of strings;
    // - Only modelID, not the index, is defined for custom/user-defined
    //   physics models;
    // - The modelID of any custom/user-defined physics model is not invariant,
    //   i.e. it is only guaranteed to be unique in any given application, but
    //   its value can be different depending on the Geant4 version, physics list
    //   and application.

    static G4bool SanityCheck();
    // This method returns "false" if "isInitialized" is false, else it does
    // a number of sanity checks, and if these are all fine, then it returns
    // "true": if any of these sanity checks are not valid, then it throws
    // a "FatalException".
    // The sanity checks performed by this method are the following:
    // it checks that the two vectors (of model IDs and model names for Geant4
    // native physics models) have the same size, the model IDs have the
    // expected values (i.e. within the allowed interval), and there are no
    // duplication of either model IDs  or model names.
    // Note that the uniqueness of the model names is checked by considering all
    // physics models - Geant4 native and custom/user-defined; there is no check
    // of the modelIDs of custom/user-defined physics models, because these
    // values are automatically assigned.
  
  private:
    G4PhysicsModelCatalog() = default;

    inline static void InsertModel( G4int modelID, G4String modelName );
  
    static G4bool isInitialized;
    static const G4int theMinAllowedModelIDValue = 10000;
    static const G4int theMaxAllowedModelIDValue = 39999;
  
    static std::vector< G4int >*    theVectorOfModelIDs;          // Non-contiguous large, positive integers
    static std::vector< G4String >* theVectorOfModelNames;        // for Geant4 native physics models
    static std::vector< G4String >* theVectorOfCustomModelNames;  // For custom/user-defined physics models  
};


inline G4int G4PhysicsModelCatalog::GetMinAllowedModelIDValue() {
  return theMinAllowedModelIDValue;
}

inline G4int G4PhysicsModelCatalog::GetMaxAllowedModelIDValue() {
  return theMaxAllowedModelIDValue;
}

inline void G4PhysicsModelCatalog::InsertModel( G4int modelID, G4String modelName ) {
  theVectorOfModelIDs->push_back( modelID );
  theVectorOfModelNames->push_back( modelName );
}

#endif
