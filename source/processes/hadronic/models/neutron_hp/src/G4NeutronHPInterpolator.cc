//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPInterpolator.hh"

  G4double G4NeutronHPInterpolator::
  GetBinIntegral(const G4InterpolationScheme & aScheme, 
                const G4double x1,const G4double x2,const G4double y1,const G4double y2)
  { // inline again later on @@@@
    G4double result = 0;
    if(aScheme==HISTO||aScheme==CHISTO||aScheme==UHISTO)
    {
      result = y1*(x2-x1);
    }
    else if(aScheme==LINLIN||aScheme==CLINLIN||aScheme==ULINLIN)
    {
      result = 0.5*(y2+y1)*(x2-x1);
    }
    else if(aScheme==LINLOG||aScheme==CLINLOG||aScheme==ULINLOG)
    {
      G4double b = (y2-y1)/(log(x2)-log(x1));
      G4double a = y1 - b*log(x1);
      result = (a-b)*(x2-x1) + b*(x2*log(x2)-x1*log(x1));
    }
    else if(aScheme==LOGLIN||aScheme==CLOGLIN||aScheme==ULOGLIN)
    {
      G4double b = (log(y2)-log(y1))/(x2-x1);
      G4double a = log(y1) - b*x1;
      result = (exp(a)/b)*(exp(b*x2)-exp(b*x1));
    }
    else if(aScheme==LOGLOG||aScheme==CLOGLOG||aScheme==ULOGLOG)
    {
      G4double b = (log(y2)-log(y1))/(log(x2)-log(x1));
      G4double a = log(y1) - b*log(x1);;
      result = (exp(a)/(b+1))*(pow(x2,b+1)-pow(x1,b+1));
    }
    else
    {
      G4Exception("Unknown interpolation scheme in G4NeutronHPVector::Integrate");
    }
    return result;
  }
  G4double G4NeutronHPInterpolator::
  GetWeightedBinIntegral(const G4InterpolationScheme & aScheme, 
                         const G4double x1,const G4double x2,const G4double y1,const G4double y2)
  { // inline again later on @@@@
    G4double result = 0;
    if(aScheme==HISTO||aScheme==CHISTO||aScheme==UHISTO)
    {
      result = 0.5*y1*(x2*x2-x1*x1);
    }
    else if(aScheme==LINLIN||aScheme==CLINLIN||aScheme==ULINLIN)
    {
      G4double b = (y2-y1)/(x2-x1);
      G4double a = y1 - b*x1;
      result = 0.5*a*(x2*x2-x1*x1) + (b/3.)*(x2*x2*x2-x1*x1*x1);
    }
    else if(aScheme==LINLOG||aScheme==CLINLOG||aScheme==ULINLOG)
    {
      G4double b = (y2-y1)/(log(x2)-log(x1));
      G4double a = y1 - b*log(x1);
      result = ( x2*x2/2. * (a-b/2.+b*log(x2)) )
              -( x1*x1/2. * (a-b/2.+b*log(x1)) );
    }
    else if(aScheme==LOGLIN||aScheme==CLOGLIN||aScheme==ULOGLIN)
    {
      G4double b = (log(y2)-log(y1))/(x2-x1);
      G4double a = log(y1) - b*x1;
      result = exp(a)/(b*b)*( exp(b*x2)*(b*x2-1.) - exp(b*x1)*(b*x1-1.) );
    }
    else if(aScheme==LOGLOG||aScheme==CLOGLOG||aScheme==ULOGLOG)
    {
      G4double b = (log(y2)-log(y1))/(log(x2)-log(x1));
      G4double a = log(y1) - b*log(x1);;
      result = exp(a)/(b+2.)*( pow(x2, b+2.) - pow(x1, b+2) );
    }
    else
    {
      G4Exception("Unknown interpolation scheme in G4NeutronHPVector::Integrate");
    }
    return result;
  }
