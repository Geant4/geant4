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
// File name:     RadmonTokenizer.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTokenizer.hh,v 1.2 2006-06-28 13:52:12 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
