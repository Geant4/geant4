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
// File name:     RadmonDetectorMultilayersLayoutCollection.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayersLayoutCollection.hh,v 1.4 2006/06/29 16:11:42 gunter Exp $
// Tag:           $Name: geant4-09-01 $
//
// Description:   Internal class to collect multilayers
//

#ifndef   RADMONDETECTORMULTILAYERSLAYOUTCOLLECTION_HH
 #define  RADMONDETECTORMULTILAYERSLAYOUTCOLLECTION_HH

 // Include files
 #include "RadmonDetectorMultilayerLayout.hh"
 #include "RadmonTLabelledCollection.hh"

 class RadmonDetectorMultilayersLayoutCollection
 {
  public:
   inline                                       RadmonDetectorMultilayersLayoutCollection();
   inline                                      ~RadmonDetectorMultilayersLayoutCollection();

   G4int                                        GetNMultilayers(void) const;
   G4bool                                       Empty(void) const;

   const RadmonDetectorMultilayerLayout &       GetMultilayer(G4int index) const;
   RadmonDetectorMultilayerLayout &             GetMultilayer(G4int index);

   G4bool                                       ExistsMultilayerByLabel(const G4String & label) const;
   G4int                                        MultiplicityMultilayerByLabel(const G4String & label) const;

   const RadmonDetectorMultilayerLayout &       FindMultilayerByLabel(const G4String & label, G4int count = 0) const;
   RadmonDetectorMultilayerLayout &             FindMultilayerByLabel(const G4String & label, G4int count = 0);

   RadmonDetectorMultilayerLayout &             CreateMultilayer(void);

   void                                         RemoveMultilayerByLabel(const G4String & label, G4int count = 0);
   void                                         RemoveMultilayersByLabel(const G4String & label);
   void                                         RemoveMultilayer(G4int index);
   void                                         RemoveAllMultilayers(void);
 
   void                                         DumpLayout(std::ostream & out, const G4String &indent=G4String()) const;

  private:
  // Hidden constructors and operators
                                                RadmonDetectorMultilayersLayoutCollection(const RadmonDetectorMultilayersLayoutCollection & copy);
   RadmonDetectorMultilayersLayoutCollection &  operator=(const RadmonDetectorMultilayersLayoutCollection & copy);

  // Private attributes
   RadmonTLabelledCollection<RadmonDetectorMultilayerLayout> multilayersCollection;
 };

 // Inline implementations
 #include "RadmonDetectorMultilayersLayoutCollection.icc"
#endif /* RADMONDETECTORMULTILAYERSLAYOUTCOLLECTION_HH */
