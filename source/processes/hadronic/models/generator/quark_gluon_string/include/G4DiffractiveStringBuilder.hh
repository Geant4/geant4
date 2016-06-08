#ifndef G4DiffractiveStringBuilder_h
#define G4DiffractiveStringBuilder_h 1

#include "globals.hh"
#include "G4KineticTrackVector.hh"
#include "G4ExcitedStringVector.hh"
#include "G4PartonPair.hh"

class G4DiffractiveStringBuilder 
    {
public:
      G4DiffractiveStringBuilder();
      G4DiffractiveStringBuilder(const G4DiffractiveStringBuilder &right);
      ~G4DiffractiveStringBuilder();
      
      int operator==(const G4DiffractiveStringBuilder &right) const;
      int operator!=(const G4DiffractiveStringBuilder &right) const;
      
      G4ExcitedString* BuildString(G4PartonPair* aParton);
      
private:     
     };
     
inline int G4DiffractiveStringBuilder::operator==(const G4DiffractiveStringBuilder &right) const
    {
    return 1;
    }

inline int G4DiffractiveStringBuilder::operator!=(const G4DiffractiveStringBuilder &right) const
    {
    return  0;
    }

#endif     
