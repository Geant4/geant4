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
// File name:     RadmonGeneratorSourceAlgorithmLayout.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorSourceAlgorithmLayout.hh,v 1.2 2006-06-28 13:53:05 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Container for the data of an algorithm of a particles source
//

#ifndef   RADMONGENERATORSOURCEALGORITHMLAYOUT_HH
 #define  RADMONGENERATORSOURCEALGORITHMLAYOUT_HH
 
 // Include files
 #include "RadmonLayoutEntityWithAttributes.hh"

 class RadmonGeneratorSourceAlgorithmLayout : public RadmonLayoutEntityWithAttributes
 {
  public:
   inline                                       RadmonGeneratorSourceAlgorithmLayout();
                                                RadmonGeneratorSourceAlgorithmLayout(const RadmonGeneratorSourceAlgorithmLayout & copy);
   inline                                      ~RadmonGeneratorSourceAlgorithmLayout();

   RadmonGeneratorSourceAlgorithmLayout &       operator=(const RadmonGeneratorSourceAlgorithmLayout & copy);

   inline const G4String &                      GetLabel(void) const;
   inline void                                  SetLabel(const G4String & label);

   inline const G4String &                      GetType(void) const;
   inline void                                  SetType(const G4String & type);

   void                                         DumpLayout(std::ostream & out, const G4String & indent = G4String()) const;

  // Private attributes
  private:
   G4String                                     algorithmLabel;
   G4String                                     algorithmType;
 };
 
 // Inline implementations
 #include "RadmonGeneratorSourceAlgorithmLayout.icc"
#endif /* RADMONGENERATORSOURCEALGORITHMLAYOUT_HH */
