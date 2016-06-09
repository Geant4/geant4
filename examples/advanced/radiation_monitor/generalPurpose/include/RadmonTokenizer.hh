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
// File name:     RadmonTokenizer.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTokenizer.hh,v 1.3 2006/06/29 16:14:33 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//
// Description:   Reads single tokens of a string
//

#ifndef   RADMONTOKENIZER_HH
 #define  RADMONTOKENIZER_HH
 
 // Include files
 #include "globals.hh"
 #include "G4String.hh"

 class RadmonTokenizer
 {
  public:
   inline                                       RadmonTokenizer(const G4String & str);
   inline                                      ~RadmonTokenizer();
   
   inline void                                  rewind(void);
   G4String                                     operator()(const char * separators=" \t\n", const char * quotations="'\"");

   inline G4bool                                eos(void) const;

  private:
  // Hidden constructors and operators
                                                RadmonTokenizer();
                                                RadmonTokenizer(const RadmonTokenizer & copy);
   RadmonTokenizer &                            operator=(const RadmonTokenizer & copy);
   
  // Private attributes
   str_size                                     position;
   G4String                                     data;
 };

 // Inline implementations 
 #include "RadmonTokenizer.icc"
#endif /* RADMONTOKENIZER_HH */
