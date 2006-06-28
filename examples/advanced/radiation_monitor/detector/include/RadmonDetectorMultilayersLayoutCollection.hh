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
// File name:     RadmonDetectorMultilayersLayoutCollection.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayersLayoutCollection.hh,v 1.3 2006-06-28 13:48:42 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
