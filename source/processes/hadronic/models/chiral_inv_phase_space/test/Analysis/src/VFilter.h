#ifndef VANAFilter_h
#define VANAFilter_h

template <class Type>
class TVANAFilter
{
  public:
  
  TVANAFilter(G4String aName) : theName(aName) {}
  virtual G4bool Accept(Type & anInput) = 0;
  virtual G4double RelativeGeometricalAcceptance() { return 1;}
  G4String GetName() {return theName;}
  
  private:
  
  TVANAFilter() {}
  
  private:
  
  G4String theName;
};

#endif
