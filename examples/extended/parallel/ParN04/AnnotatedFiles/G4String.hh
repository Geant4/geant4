/*
  This file only serves as a place-holder for the annotations of G4String class,
  which is a type alias in Geant4. 
  
  Only the marshaling code (MarshaledG4String.hh) generated from this file
   should be used in your code. Please make sure that this G4String.hh
   itself is NOT included in your code. Otherwise, it will override Geant4's
   valid G4String class and bad things would happen.

  vietha 2003.05.08
*/

#ifndef G4String_h
#define G4String_h

//MSH_BEGIN
class G4String{
    int dummy; /*MSH: manual
    { memcpy($$, param->c_str(), param->size());
    *($$+param->size()) = '\0'; 
    }
    { G4String* s = new G4String($$);
       memcpy(param, s, sizeof(G4String));
     }
     { int size = param->size()+1;
       while(size%8) size++;
       $SIZE = size; } */
};
//MSH_END
#endif
