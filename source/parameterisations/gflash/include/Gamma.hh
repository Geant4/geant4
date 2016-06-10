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
//
// $Id: Gamma.hh 68057 2013-03-13 14:46:00Z gcosmo $
//
//
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  Gamma
//
//  Class description:
//
//  Gamma function.

//---------------------------------------------------------------
#ifndef MYGAMMA_H
#define MYGAMMA_H

#include <cmath>

class MyGamma
{
  public:

    MyGamma ();
    ~MyGamma();
  
    double Gamma(double z);
    double Gamma(double a,double x);

  private:

    double GamCf(double a,double x);
    double GamSer(double a,double x);
  
    // Abs
    static short  Abs(short d)  { return (d > 0) ? d : -d; }
    static int    Abs(int d)    { return (d > 0) ? d : -d; }
    static long   Abs(long d)   { return (d > 0) ? d : -d; }
    static float  Abs(float d)  { return (d > 0) ? d : -d; }
    static double Abs(double d) { return (d > 0) ? d : -d; }
    static double LnGamma(double z); 
    static double Log(double x) { return std::log(x); }
    static double Exp(double x) { return std::exp(x); }
};

#endif
