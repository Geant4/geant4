#ifndef G4SoftStringBuilder_h
#define G4SoftStringBuilder_h 1

#include "globals.hh"
#include "G4KineticTrackVector.hh"
#include "G4ExcitedStringVector.hh"
#include "G4PartonPair.hh"


class G4SoftStringBuilder 
    {
public:
    G4SoftStringBuilder();
    G4SoftStringBuilder(const G4SoftStringBuilder &right);
    ~G4SoftStringBuilder();

    int operator==(const G4SoftStringBuilder &right) const;
    int operator!=(const G4SoftStringBuilder &right) const;

    G4ExcitedString* BuildString(G4PartonPair * aPair);        

private:     
    };
     
inline int G4SoftStringBuilder::operator==(const G4SoftStringBuilder &right) const
    {
    return 1;
    }

inline int G4SoftStringBuilder::operator!=(const G4SoftStringBuilder &right) const
    {
    return  0;
    }

#endif     
