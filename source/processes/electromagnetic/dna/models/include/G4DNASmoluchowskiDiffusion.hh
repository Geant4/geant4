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
/*
 * G4DNASmoluchowskiDiffusion.hh
 *
 *  Created on: 2 févr. 2015
 *      Author: matkara
 */

#ifndef G4DNASMOLUCHOWSKIDIFFUSION_HH_
#define G4DNASMOLUCHOWSKIDIFFUSION_HH_

//#if __cplusplus >= 201103L
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

//#define DNADEV_TEST

#ifdef DNADEV_TEST
#include <TF1.h>
#endif

#include <cassert>

#ifndef DNADEV_TEST
#include "globals.hh"
#include "Randomize.hh"
#endif

#ifdef DNADEV_TEST
#include "TRandom.h"
TRandom root_random;
double G4UniformRand()
{
  return root_random.Rndm();
}

#define G4cout std::cout
#define G4endl std::endl
#endif

class G4DNASmoluchowskiDiffusion
{
public:
  G4DNASmoluchowskiDiffusion(double epsilon = 1e-5);
  virtual ~G4DNASmoluchowskiDiffusion();

  static double ComputeS(double r, double D, double t)
  {
    double sTransform = r / (2. * std::sqrt(D * t));
    return sTransform;
  }

  static double ComputeDistance(double sTransform, double D, double t)
  {
    return sTransform * 2. * std::sqrt(D * t);
  }

  static double ComputeTime(double sTransform, double D, double r)
  {
    return std::pow(r / sTransform, 2.) / (4. * D);
  }

  //====================================================

  double GetRandomDistance(double _time, double D)
  {
    double proba = G4UniformRand();
//    G4cout << "proba = " << proba << G4endl;
    double sTransform = GetInverseProbability(proba);
//    G4cout << "sTransform = " << sTransform << G4endl;
    return ComputeDistance(sTransform, D, _time);
  }

  double GetRandomTime(double distance, double D)
  {
    double proba = G4UniformRand();
    double sTransform = GetInverseProbability(proba);
    return ComputeTime(sTransform, D, distance);
  }

  double EstimateCrossingTime(double proba,
                              double distance,
                              double D)
  {
    double sTransform = GetInverseProbability(proba);
    return ComputeTime(sTransform, D, distance);
  }

  //====================================================
  // 1-value transformation

  // WARNING : this is NOT the differential probability
  // this is the derivative of the function GetCumulativeProbability
  static double GetDifferential(double sTransform)
  {
    static double constant = -4./std::sqrt(3.141592653589793);
    return sTransform*sTransform*std::exp(-sTransform*sTransform)*constant; // -4*sTransform*sTransform*exp(-sTransform*sTransform)/sqrt(3.141592653589793)
  }

  static double GetDensityProbability(double r, double _time, double D)
  {
    static double my_pi = 3.141592653589793;
    static double constant = 4.*my_pi/std::pow(4.*my_pi, 1.5);
    return r*r/std::pow(D * _time,1.5)*std::exp(-r*r/(4. * D * _time))*constant;
  }

  //====================================================
  // BOUNDING BOX
  struct BoundingBox
  {
    double fXmax;
    double fXmin;
    double fXmaxDef;
    double fXminDef;
    double fToleranceY;
    double fSum;
    double    fIncreasingCumulativeFunction;

    enum PreviousAction
    {
      IncreaseProba,
      DecreaseProba,
      Undefined
    };

    PreviousAction fPreviousAction;

    BoundingBox(double xmin,
                double xmax,
                double toleranceY) :
     fXmax(xmax), fXmin(xmin),
     fToleranceY(toleranceY),
     fSum(0)
    {
      if(fXmax < fXmin)
      {
        double tmp = fXmin;
        fXmin = fXmax;
        fXmax = tmp;
      }
      
      fXminDef = fXmin;
      fXmaxDef = fXmax;
      fPreviousAction = BoundingBox::Undefined;
      fIncreasingCumulativeFunction = (GetCumulativeProbability(fXmax) - GetCumulativeProbability(fXmin))/(fXmax-fXmin);
    }
    
    void Print()
    {
      G4cout << "fXmin: " << fXmin << " | fXmax: " << fXmax << G4endl;
    }

    bool Propose(double proposedXValue,
                 double proposedProba,
                 double nextProba,
                 double& returnedValue)
    {
//      G4cout << "---------------------------" << G4endl;
//      G4cout << "Proposed x value: " << proposedXValue
//          << "| proposedProba: " << proposedProba
//          << "| nextProba: " << nextProba
//          << " | fXmin: " << fXmin << " (" << G4DNASmoluchowskiDiffusion::GetCumulativeProbability(fXmin) <<")"
//          << " | fXmax: " << fXmax << " (" << G4DNASmoluchowskiDiffusion::GetCumulativeProbability(fXmax) <<")"
//          << G4endl;

      bool returnFlag = false;
      
      if(proposedProba < nextProba-fToleranceY) // proba trop petite ==> augmente
      {
        // G4cout << "proposedProba < nextProba-fToleranceY" << G4endl;

        if(fIncreasingCumulativeFunction > 0) // croissant
        {
          if(proposedXValue > fXmin)
            fXmin = proposedXValue;
        }
        else if(fIncreasingCumulativeFunction < 0) // decroissant
        {
          if(proposedXValue < fXmax)
            fXmax = proposedXValue;
        }
        
        returnedValue = (fXmax + fXmin)/2;
        returnFlag = false;
        fPreviousAction = BoundingBox::IncreaseProba;
      }
      else if(proposedProba > nextProba+fToleranceY) // proba trop grande
      {
        // G4cout << "proposedProba > nextProba+fToleranceY" << G4endl;

        if(fIncreasingCumulativeFunction>0)
        {
          if(proposedXValue < fXmax)
            fXmax = proposedXValue;
        }
        else if(fIncreasingCumulativeFunction<0)
        {
          if(proposedXValue > fXmin)
          {
            fXmin = proposedXValue;
          }
        }
        
        returnedValue = (fXmax + fXmin)/2;
        returnFlag = false;
        fPreviousAction = BoundingBox::DecreaseProba;
      }
      else
      {
        // G4cout << "IN THE INTERVAL !! : " << nextProba << G4endl;
        fSum = proposedProba;
        
        // Assuming search for next proba is increasing
        if(fIncreasingCumulativeFunction<0)
        {
         fXmin = fXminDef;
         fXmax = proposedXValue;
        }
        else if(fIncreasingCumulativeFunction>0)
        {
          fXmin = proposedXValue;
          fXmax = fXmaxDef;
        }
        returnFlag = true;
        fPreviousAction = BoundingBox::Undefined;
      }
      
      return returnFlag;
    }
  };
  // END OF BOUNDING BOX
  //==============================
  
  void PrepareReverseTable(double xmin, double xmax)
  {
    double x = xmax;
    int index = 0;
    double nextProba = fEpsilon;
    double proposedX;

    BoundingBox boundingBox(xmin, xmax, fEpsilon*1e-5);

    while(index <= fNbins)
    // in case GetCumulativeProbability is exact (digitally speaking), replace with:
    // while(index <= fNbins+1)
    {
      nextProba = fEpsilon*index;

      double newProba = GetCumulativeProbability(x);

      if(boundingBox.Propose(x, newProba, nextProba, proposedX))
      {
        fInverse[index] = x;
        index++;
      }
      else
      {
        if(x == proposedX)
        {
          G4cout << "BREAK : x= " << x << G4endl;
          G4cout << "index= " << index << G4endl;
          G4cout << "nextProba= " << nextProba << G4endl;
          G4cout << "newProba= " << newProba << G4endl;
          abort();
        }
        x = proposedX;
      }
    }
    
    fInverse[fNbins+1] = 0; // P(1) = 0, because we want it exact !
    // Tips to improve the exactness: get an better value of pi, get better approximation of erf and exp, use long double instead of double
//    boundingBox.Print();
  }

  static double GetCumulativeProbability(double sTransform)
  {
    static double constant = 2./std::sqrt(3.141592653589793);
    return erfc(sTransform) + constant*sTransform*std::exp(-sTransform*sTransform);
  }

  double GetInverseProbability(double proba) // returns sTransform
  {
    size_t index_low = (size_t) trunc(proba/fEpsilon);
    
    if(index_low == 0) // assymptote en 0
    {
      index_low = 1;
      size_t index_up = 2;
      double low_y = fInverse[index_low];
      double up_y = fInverse[index_up];
      double low_x = index_low*fEpsilon;
      double up_x = proba+fEpsilon;
      double tangente = (low_y-up_y)/(low_x - up_x); // ou utiliser GetDifferential(proba) ?
      // double tangente = GetDifferential(proba);
      return low_y + tangente*(proba-low_x);
    }

    size_t index_up = index_low+1;
    if(index_low > fInverse.size()) return fInverse.back();
    double low_y = fInverse[index_low];
    double up_y = fInverse[index_up];

    double low_x = index_low*fEpsilon;
    double up_x = low_x+fEpsilon;

    if(up_x > 1) // P(1) = 0
    {
      up_x = 1;
      up_y = 0; // more general : fInverse.back()
    }

    double tangente = (low_y-up_y)/(low_x - up_x);

    return low_y + tangente*(proba-low_x);
  }

  double PlotInverse(double* x, double* )
  {
    return GetInverseProbability(x[0]);
  }

  double Plot(double* x, double* )
  {
    return GetDifferential(x[0]);
  }


  void InitialiseInverseProbability(double xmax = 3e28)
  {
    // x > x'
    // P'(x) = p(x') = lim(x->x') (P(x') - P(x))/(x'-x)
    // x'-x = (P(x') - P(x))/p(x')
    // x = x' - (P(x') - P(x))/p(x')

    // fInverse initialized in the constructor

    assert(fNbins !=0);
    PrepareReverseTable(0,xmax);
  }

  std::vector<double> fInverse;
  int fNbins;
  double fEpsilon;
};

#endif /* SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MODELS_G4DNASMOLUCHOWSKIDIFFUSION_HH_ */
