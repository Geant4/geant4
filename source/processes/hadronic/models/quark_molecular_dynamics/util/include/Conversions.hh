#ifndef __CONVERSIONS__
#define __CONVERSIONS__

#include "Error.hh"

template <class t>
t& Conversion(char* string,t& value);

class ConversionImpossible : public Error
{
public:
   ConversionImpossible(char* Type,char* string);
   ~ConversionImpossible();
   void writeMessage(G4std::ostream&) const;
private:
   char * String,* Type;
};

#endif
