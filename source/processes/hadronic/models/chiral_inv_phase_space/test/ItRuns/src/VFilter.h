#ifndef VANAFilter_h
#define VANAFilter_h

template <class Type>
class TVANAFilter
{
  public:
  
  virtual G4bool Accept(Type & anInput) = 0;
};

#endif
