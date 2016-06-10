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
// $Id: G4SimplexDownhill.hh 79763 2014-03-13 15:16:47Z gcosmo $
//
// Class description:
//
// Class implementing minimization of a function of n variables.
// Reference: "A Simplex method for function minimization"
//            by J. A. Nelder and R. Mead, Computer Journal, 7, 308 (1965)
// and also: "Numerical Recipes in C: the art of scientific computing"
//            by William H., Cambridge University Press ISBN 0521437202 (1992)

// Author: Tatsumi Koi (SLAC/SCCS), 2007
// --------------------------------------------------------------------------

#ifndef G4SimplexDownhill_hh
#define G4SimplexDownhill_hh

#include "globals.hh"

#include <vector> 
#include <algorithm> 

template<class T>
class G4SimplexDownhill
{

   public: // with description

      G4SimplexDownhill( T* tp , G4int n )
        : currentValue(0.), target(tp), numberOfVariable(n)
      { init(); }

      ~G4SimplexDownhill();

      G4double GetMinimum();

      std::vector< G4double > GetMinimumPoint();


   private:

      G4double getValue( std::vector< G4double > x )
      { return target->GetValueOfMinimizingFunction( x ); }

      void initialize();
      std::vector< std::vector< G4double > > currentSimplex;

      void calHeights();
      std::vector< G4double > currentHeights;
      G4double currentValue;

      std::vector< G4double > calCentroid( G4int );

      G4bool isItGoodEnough();

      std::vector< G4double > getReflectionPoint( std::vector< G4double > ,
                                                  std::vector< G4double > );
      std::vector< G4double > getExpansionPoint( std::vector< G4double > ,
                                                 std::vector< G4double > );
      std::vector< G4double > getContractionPoint( std::vector< G4double > ,
                                                   std::vector< G4double > );

      void doDownhill();

      void init();

   private:

      T* target;

      G4int numberOfVariable; 

      G4double alpha;
      G4double beta;
      G4double gamma;
      G4double max_se;
      G4double max_ratio;
      G4int maximum_no_trial;
      G4bool minimized;

      std::vector< G4double > minimumPoint;
};

#include "G4SimplexDownhill.icc"

#endif
