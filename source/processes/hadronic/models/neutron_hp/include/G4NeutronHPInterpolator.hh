// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPInterpolator.hh,v 1.5 1999-12-15 14:53:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPInterpolator_h
#define G4NeutronHPInterpolator_h 1

#include "globals.hh"
#include "G4InterpolationScheme.hh"
#include "Randomize.hh"
#include "G4ios.hh"


class G4NeutronHPInterpolator
{
  public:
  
  G4NeutronHPInterpolator(){}
  ~G4NeutronHPInterpolator()
   {
    //  G4cout <<"deleted the interpolator"<<G4endl;
   }
  
  inline G4double Lin(G4double x,G4double x1,G4double x2,G4double y1,G4double y2)
  {
    G4double slope=0, off=0;
    slope = (y2-y1)/(x2-x1);
    off = y2-x2*slope;
    G4double y = x*slope+off;
    return y;
  }
  
  inline G4double Interpolate(G4InterpolationScheme aScheme,
                              G4double x, G4double x1, G4double x2,
                                          G4double y1, G4double y2) const;       
  
  G4double 
  GetBinIntegral(const G4InterpolationScheme & aScheme, 
                const G4double x1,const G4double x2,const G4double y1,const G4double y2);

  G4double 
  GetWeightedBinIntegral(const G4InterpolationScheme & aScheme, 
                const G4double x1,const G4double x2,const G4double y1,const G4double y2);

  private:
  
  inline G4double Histogram(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const;
  inline G4double LinearLinear(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const;
  inline G4double LinearLogarithmic(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const;
  inline G4double LogarithmicLinear(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const;
  inline G4double LogarithmicLogarithmic(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const;
  inline G4double Random(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const;

};

inline G4double G4NeutronHPInterpolator::
Interpolate(G4InterpolationScheme aScheme,
            G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  G4double result;
  G4int theScheme = aScheme;
  theScheme = theScheme%CSTART;
  switch(theScheme)
  {
    case 1:
      result = Histogram(x, x1, x2, y1, y2);
      break;
    case 2:
      result = LinearLinear(x, x1, x2, y1, y2);
      break;
    case 3:
      result = LinearLogarithmic(x, x1, x2, y1, y2);
      break;
    case 4:
      result = LogarithmicLinear(x, x1, x2, y1, y2);
      break;
    case 5:
      result = LogarithmicLogarithmic(x, x1, x2, y1, y2);
      break;
    case 6:
      result = Random(x, x1, x2, y1, y2);
      break;
    default:
      G4cout << "theScheme = "<<theScheme<<G4endl;
      G4Exception("G4NeutronHPInterpolator::Carthesian Invalid InterpolationScheme");
      break;
  }
  return result;
}

inline G4double G4NeutronHPInterpolator::
Histogram(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  G4double result;
  result = y1;
  return result;
}

inline G4double G4NeutronHPInterpolator::
LinearLinear(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  G4double slope=0, off=0;
  slope = (y2-y1)/(x2-x1);
  off = y2-x2*slope;
  G4double y = x*slope+off;
  return y;
}

inline G4double G4NeutronHPInterpolator::
LinearLogarithmic(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  G4double result;
  result = LinearLinear(log(x), log(x1), log(x2), y1, y2);
  return result;
}
  
inline G4double G4NeutronHPInterpolator::
LogarithmicLinear(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  G4double result;
  result = LinearLinear(x, x1, x2, log(y1), log(y2));
  result = exp(result);
  return result;
}
  
inline G4double G4NeutronHPInterpolator::
LogarithmicLogarithmic(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  G4double result;
  result = LinearLinear(log(x), log(x1), log(x2), log(y1), log(y2));
  result = exp(result);
  return result;
}

inline G4double G4NeutronHPInterpolator::
Random(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  G4double result;
  result = y1+G4UniformRand()*(y2-y1);
  return result;
}

#endif
