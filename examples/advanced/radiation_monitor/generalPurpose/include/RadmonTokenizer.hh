//
// File name:     RadmonTokenizer.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTokenizer.hh,v 1.1 2005-09-14 12:30:15 capra Exp $
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
