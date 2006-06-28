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
// File name:     RadmonDetectorEnvironmentLayout.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorEnvironmentLayout.hh,v 1.4 2006-06-28 13:46:42 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Internal class to describe the environment layout
//

#ifndef   RADMONDETECTORENVIRONMENTLAYOUT_HH
 #define  RADMONDETECTORENVIRONMENTLAYOUT_HH
 
 // Include files
 #include "RadmonLayoutEntityWithAttributes.hh"
 
 class RadmonDetectorEnvironmentLayout : public RadmonLayoutEntityWithAttributes
 {
  public:
   inline                                     RadmonDetectorEnvironmentLayout();
   inline                                    ~RadmonDetectorEnvironmentLayout();

   inline G4bool                              IsEnabled(void) const;
   inline void                                Enable(void);
   inline void                                Disable(void);

   inline const G4String &                    GetType(void) const;
   inline void                                SetType(const G4String & type);

   void                                       DumpLayout(std::ostream & out, const G4String & indent=G4String("")) const;

  private:
                                              RadmonDetectorEnvironmentLayout(const RadmonDetectorEnvironmentLayout & copy);
   RadmonDetectorEnvironmentLayout &          operator=(const RadmonDetectorEnvironmentLayout & copy);

   G4bool                                     enabled;
   G4String                                   environmentType;
 };

 // Inline implementations
 #include "RadmonDetectorEnvironmentLayout.icc"
#endif /* RADMONDETECTORENVIRONMENTLAYOUT_HH */

