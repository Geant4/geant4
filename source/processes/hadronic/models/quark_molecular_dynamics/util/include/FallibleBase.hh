#ifndef _FallibleBase_H
#define _FallibleBase_H

#include "Definitions.hh"
#include "Boolean.h"
#include "Error.hh"

class G4std::ostream;

class FallibleBase {
  public:
    FallibleBase(Boolean state) : valid(state) {}
    Boolean isInvalid() const { return !valid; }
    Boolean isValid() const { return valid; }
    void invalidate() { valid = Boolean::False; }
    void validate() { valid = Boolean::True; }
    
    class ErrUsedInInvalidState
      : public Error
    {
      public:
        virtual void writeMessage(G4std::ostream& os) const;
    };

  protected:
    void throwError() const { Throw(ErrUsedInInvalidState()); }
  private:
    Boolean valid;
};

#endif // _FallibleBase_H


