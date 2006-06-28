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
// File name:     RadmonDetectorLayerLayout.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerLayout.hh,v 1.4 2006-06-28 13:47:33 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Internal class that describes layer informations
//

#ifndef   RADMONDETECTORLAYERLAYOUT_HH
 #define  RADMONDETECTORLAYERLAYOUT_HH

 // Include files
 #include "RadmonLayoutEntityWithAttributes.hh"

 class RadmonDetectorLayerLayout : public RadmonLayoutEntityWithAttributes
 {
  public:
   inline                                       RadmonDetectorLayerLayout();
   inline                                       RadmonDetectorLayerLayout(const RadmonDetectorLayerLayout & copy);
   inline                                      ~RadmonDetectorLayerLayout();

   inline RadmonDetectorLayerLayout &           operator=(const RadmonDetectorLayerLayout & copy);

   inline const G4String &                      GetLabel(void) const;
   inline const G4String &                      GetType(void) const;
   inline G4double                              GetThickness(void) const;
   inline void                                  SetLabel(const G4String & label);
   inline void                                  SetType(const G4String & type);
   inline void                                  SetThickness(G4double thickness);

   void                                         DumpLayout(std::ostream & out, const G4String &indent=G4String()) const;

  private:
  // Private attributes
   G4String                                     layerLabel;
   G4String                                     layerType;
   G4double                                     layerThickness;
 };

 // Inline implementations
 #include "RadmonDetectorLayerLayout.icc"
#endif /* RADMONDETECTORLAYERLAYOUT_HH */
