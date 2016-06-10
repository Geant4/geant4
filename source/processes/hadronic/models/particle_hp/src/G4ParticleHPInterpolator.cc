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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPInterpolator.hh"
#include "G4Pow.hh"

  G4double G4ParticleHPInterpolator::
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
      if(x1==0) result = y1;
      else if(x2==0) result = y2;
      else
      {
        G4double b = (y2-y1)/(G4Log(x2)-G4Log(x1));
        G4double a = y1 - b*G4Log(x1);
        result = (a-b)*(x2-x1) + b*(x2*G4Log(x2)-x1*G4Log(x1));
      }
    }
    else if(aScheme==LOGLIN||aScheme==CLOGLIN||aScheme==ULOGLIN)
    {
      if ( y1==0||y2==0 ) { 
         result =0;
      } else {
        // G4double b = (std::log(y2)-std::log(y1))/(x2-x1);
        // G4double a = std::log(y1) - b*x1;
        //***************************************************************
        //EMendoza:
        //result = (std::exp(a)/b)*(std::exp(b*x2)-std::exp(b*x1));
        //***************************************************************
        if ( y1!=y2 ) {
           result = (x2-x1)*(y2-y1)/G4Log(y2/y1);
        } else { 
          result = y2*(x2-x1);
        }
        //***************************************************************
      }
    }
    else if(aScheme==LOGLOG||aScheme==CLOGLOG||aScheme==ULOGLOG)
    {
      if(x1==0) result = y1;
      else if(x2==0) result = y2;
      else if(y1==0||y2==0) result =0;
      else
      {      
        G4double b = (G4Log(y2)-G4Log(y1))/(G4Log(x2)-G4Log(x1));
        G4double a = G4Log(y1) - b*G4Log(x1);;
        result = (G4Exp(a)/(b+1))*(G4Pow::GetInstance()->powA(x2,b+1)-G4Pow::GetInstance()->powA(x1,b+1));
      }
    }
    else
    {
      throw G4HadronicException(__FILE__, __LINE__, "Unknown interpolation scheme in G4ParticleHPVector::Integrate");
    }
    return result;
  }
  G4double G4ParticleHPInterpolator::
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
      //        G4double b = (y2-y1)/(x2-x1);
      //        G4double a = y1 - b*x1;
      //        result = 0.5*a*(x2*x2-x1*x1) + (b/3.)*(x2*x2*x2-x1*x1*x1);
      //  Factor out x2-x1 to avoid divide by zero
      
      result = (y1*x2 - y2*x1)*(x2 + x1)/2. + (y2-y1)*(x2*x2 + x2*x1 + x1*x1)/3.;
    }
    else if(aScheme==LINLOG||aScheme==CLINLOG||aScheme==ULINLOG)
    {
      if(x1==0) result = y1;
      else if(x2==0) result = y2;
      else
      {
        G4double b = (y2-y1)/(G4Log(x2)-G4Log(x1));
        G4double a = y1 - b*G4Log(x1);
        result = ( x2*x2/2. * (a-b/2.+b*G4Log(x2)) )
                -( x1*x1/2. * (a-b/2.+b*G4Log(x1)) );
      }
    }
    else if(aScheme==LOGLIN||aScheme==CLOGLIN||aScheme==ULOGLIN)
    {
      if(y1==0||y2==0) result = 0;
      else
      {
        G4double b = (G4Log(y2)-G4Log(y1))/(x2-x1);
        G4double a = G4Log(y1) - b*x1;
        result = G4Exp(a)/(b*b)*( G4Exp(b*x2)*(b*x2-1.) - G4Exp(b*x1)*(b*x1-1.) );
      }
    }
    else if(aScheme==LOGLOG||aScheme==CLOGLOG||aScheme==ULOGLOG)
    {
      if(x1==0) result = y1;
      else if(x2==0) result = y2;
      else if(y1==0||y2==0) result = 0;
      else
      {
        G4double b = (G4Log(y2)-G4Log(y1))/(G4Log(x2)-G4Log(x1));
        G4double a = G4Log(y1) - b*G4Log(x1);;
        result = G4Exp(a)/(b+2.)*( G4Pow::GetInstance()->powA(x2, b+2.) - G4Pow::GetInstance()->powA(x1, b+2) );
      }
    }
    else
    {
      throw G4HadronicException(__FILE__, __LINE__, "Unknown interpolation scheme in G4ParticleHPVector::Integrate");
    }
    return result;
  }
