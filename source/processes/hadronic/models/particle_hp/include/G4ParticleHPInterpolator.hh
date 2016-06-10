//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// 080809 Change interpolation scheme of "histogram", now using LinearLinear
//        For multidimensional interpolations By T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPInterpolator_h
#define G4ParticleHPInterpolator_h 1

#include "globals.hh"
#include "G4InterpolationScheme.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4HadronicException.hh"


class G4ParticleHPInterpolator
{
  public:
  
  G4ParticleHPInterpolator(){}
  ~G4ParticleHPInterpolator()
   {
    //  G4cout <<"deleted the interpolator"<<G4endl;
   }
  
  inline G4double Lin(G4double x,G4double x1,G4double x2,G4double y1,G4double y2)
  {
    G4double slope=0, off=0;
    if(x2-x1==0) return (y2+y1)/2.;
    slope = (y2-y1)/(x2-x1);
    off = y2-x2*slope;
    G4double y = x*slope+off;
    return y;
  }
  
  inline G4double Interpolate(G4InterpolationScheme aScheme,
                              G4double x, G4double x1, G4double x2,
                                          G4double y1, G4double y2) const;       
  inline G4double Interpolate2(G4InterpolationScheme aScheme,
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

inline G4double G4ParticleHPInterpolator::
Interpolate(G4InterpolationScheme aScheme,
            G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  G4double result(0);
  G4int theScheme = aScheme;
  theScheme = theScheme%CSTART_;
  switch(theScheme)
  {
    case 1:
      //080809
      //result = Histogram(x, x1, x2, y1, y2);
      result = LinearLinear(x, x1, x2, y1, y2);
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
      throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPInterpolator::Carthesian Invalid InterpolationScheme");
      break;
  }
  return result;
}

inline G4double G4ParticleHPInterpolator::
Interpolate2(G4InterpolationScheme aScheme,
            G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  G4double result(0);
  G4int theScheme = aScheme;
  theScheme = theScheme%CSTART_;
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
      throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPInterpolator::Carthesian Invalid InterpolationScheme");
      break;
  }
  return result;
}

inline G4double G4ParticleHPInterpolator::
Histogram(G4double , G4double , G4double , G4double y1, G4double ) const
{
  G4double result;
  result = y1;
  return result;
}

inline G4double G4ParticleHPInterpolator::
LinearLinear(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  G4double slope=0, off=0;
  if(x2-x1==0) return (y2+y1)/2.;
  slope = (y2-y1)/(x2-x1);
  off = y2-x2*slope;
  G4double y = x*slope+off;
  return y;
}

inline G4double G4ParticleHPInterpolator::
LinearLogarithmic(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  G4double result;
  if(x==0) result = y1+y2/2.;
  else if(x1==0) result = y1;
  else if(x2==0) result = y2;
  else result = LinearLinear(G4Log(x), G4Log(x1), G4Log(x2), y1, y2);
  return result;
}
  
inline G4double G4ParticleHPInterpolator::
LogarithmicLinear(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  G4double result;
  if(y1==0||y2==0) result = 0;
  else 
  {
    result = LinearLinear(x, x1, x2, G4Log(y1), G4Log(y2));
    result = G4Exp(result);
  }
  return result;
}
  
inline G4double G4ParticleHPInterpolator::
LogarithmicLogarithmic(G4double x, G4double x1, G4double x2, G4double y1, G4double y2) const
{
  if(x==0) return y1+y2/2.;
  else if(x1==0) return y1;
  else if(x2==0) return y2;
  G4double result;
  if(y1==0||y2==0) result = 0;
  else 
  {
    result = LinearLinear(G4Log(x), G4Log(x1), G4Log(x2), G4Log(y1), G4Log(y2));
    result = G4Exp(result);
  }
  return result;
}

inline G4double G4ParticleHPInterpolator::
Random(G4double , G4double , G4double , G4double y1, G4double y2) const
{
  G4double result;
  result = y1+G4UniformRand()*(y2-y1);
  return result;
}

#endif
