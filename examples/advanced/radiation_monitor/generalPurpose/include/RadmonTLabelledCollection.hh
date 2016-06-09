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
// File name:     RadmonTLabelledCollection.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTLabelledCollection.hh,v 1.4 2006/06/29 16:14:29 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//
// Description:   Container for objects with the following methods
//                 - Constructor
//                 - Copy-constructor
//                 - Destructor
//                 - Copy operation
//                 - const G4String & GetLabel(void) const
//

#ifndef   RADMONTLABELLEDCOLLECTION_HH
 #define  RADMONTLABELLEDCOLLECTION_HH
 
 // Include files
 #include "globals.hh"

 #include <vector>
 
 // Forward declarations
 class G4String;

 template <typename LabelledObject>
 class RadmonTLabelledCollection
 {
  public:
                                                RadmonTLabelledCollection();
                                                RadmonTLabelledCollection(const RadmonTLabelledCollection<LabelledObject> & copy);
                                               ~RadmonTLabelledCollection();

   RadmonTLabelledCollection<LabelledObject> &  operator=(const RadmonTLabelledCollection<LabelledObject> & copy);

   G4int                                        GetNItems() const;
   G4bool                                       Empty() const;

   G4bool                                       ExistsItemByLabel(const G4String & label) const;
   G4int                                        MultiplicityItemByLabel(const G4String & label) const;

   const LabelledObject &                       GetItem(G4int index) const;
   LabelledObject &                             GetItem(G4int index);
   const LabelledObject &                       FindItemByLabel(const G4String & label, G4int count = 0) const;
   LabelledObject &                             FindItemByLabel(const G4String & label, G4int count = 0);

   LabelledObject &                             InsertItemAfter(G4int index);
   LabelledObject &                             InsertItemBefore(G4int index);
   LabelledObject &                             AppendItem(void);
   LabelledObject &                             PrependItem(void);

   void                                         RemoveItemByLabel(const G4String & label, G4int count = 0);
   void                                         RemoveItemsByLabel(const G4String & label);
   void                                         RemoveItem(G4int index);
   void                                         RemoveItemsByRange(G4int first, G4int last);
   void                                         RemoveAllItems();

  private:
   typedef std::vector<LabelledObject>          Collection;
   Collection                                   collection;
 };

 // Inline implementations 
 #include "RadmonTLabelledCollection.icc"
#endif /* RADMONTLABELLEDCOLLECTION_HH */
