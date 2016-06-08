#ifndef G4QGSMParameters_h
#define G4QGSMParameters_h 1

#include "globals.hh"

class G4QGSMParameters 
    {
public:
      G4QGSMParameters();
      G4QGSMParameters(const G4QGSMParameters &right);
      ~G4QGSMParameters();
      
      int operator==(const G4QGSMParameters &right) const;
      int operator!=(const G4QGSMParameters &right) const;
private:     
     };
     
inline int G4QGSMParameters::operator==(const G4QGSMParameters &right) const
    {
    return 1;
    }

inline int G4QGSMParameters::operator!=(const G4QGSMParameters &right) const
    {
    return  0;
    }

#endif     
