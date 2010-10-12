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
// $Id: G4ConvergenceTester.hh,v 1.1 2010-10-12 07:24:27 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class description:
//
// Convergence Tests for Monte Carlo results.
//
// Reference
// MCNP(TM) -A General Monte Carlo N-Particle Transport Code
// Version 4B
// Judith F. Briesmeister, Editor
// LA-12625-M, Issued: March 1997, UC 705 and UC 700
// CHAPTER 2. GEOMETRY, DATA, PHYSICS, AND MATHEMATICS
//        VI. ESTIMATION OF THE MONTE CARLO PRECISION
//
// Positives numbers are assumed for input values

// Author: Tatsumi Koi (SLAC/SCCS)
// 
// --------------------------------------------------------------------

#ifndef G4ConvergenceTester
#define G4ConvergenceTester_h 1

#include "G4SimplexDownhill.hh"

#include "G4Timer.hh"
#include "globals.hh"

#include <map>
#include <vector>

class G4ConvergenceTester 
{
   public:

      G4ConvergenceTester();
     ~G4ConvergenceTester();
      G4ConvergenceTester( G4double );

   public:

      void AddScore( G4double );

      void ShowHistory();
      void ShowResult();

      inline G4double GetValueOfMinimizingFunction( std::vector<G4double> x )
             { return slope_fitting_function( x ); }

   private:

      void calStat();

      G4double GetMean() { return mean; }
      G4double GetStandardDeviation() { return sd; }
      G4double GetVariance() { return var; }
      G4double GetR() { return r; }
      G4double GetEfficiency() { return efficiency; }
      G4double GetR2eff() { return r2eff; }
      G4double GetR2int() { return r2int; }
      G4double GetShift() { return shift; }
      G4double GetVOV() { return vov; }
      G4double GetFOM() { return fom; }

      void calc_grid_point_of_history();
      void calc_stat_history();
      void check_stat_history();
      G4double calc_Pearson_r( G4int, std::vector<G4double>,
                               std::vector<G4double> );
      G4bool is_monotonically_decrease( std::vector<G4double> ); 
      void calc_slope_fit( std::vector< G4double > );
      G4double slope_fitting_function( std::vector< G4double > );

   private:

      std::map< G4int , G4double > nonzero_histories;
        // (ith-history , score value)
      G4int n;
        // number of history
      G4double sum; // sum of scores;

      G4Timer* timer; 
      std::vector<G4double> cpu_time; 

      G4double mean; 
      G4double var; 
      G4double sd; 
      G4double r;          // relative err sd/mean/sqrt(n) 
      G4double efficiency; // rate of non zero score 
      G4double r2eff; 
      G4double r2int; 
      G4double shift; 
      G4double vov; 
      G4double fom; 

      G4double largest;
      G4int largest_score_happened;
 
      G4double mean_1; 
      G4double var_1; 
      G4double sd_1; 
      G4double r_1;        // relative err sd/mean/sqrt(n) 
      G4double shift_1; 
      G4double vov_1; 
      G4double fom_1; 

      G4int noBinOfHistory;
      std::vector< G4int > history_grid;
      std::vector< G4double > mean_history;
      std::vector< G4double > var_history;
      std::vector< G4double > sd_history;
      std::vector< G4double > r_history;
      std::vector< G4double > vov_history;
      std::vector< G4double > fom_history;
      std::vector< G4double > shift_history;
      std::vector< G4double > e_history;
      std::vector< G4double > r2eff_history;
      std::vector< G4double > r2int_history;

      G4double slope; 
      std::vector< G4double > largest_scores; 
      std::vector< G4double > f_xi;
      std::vector< G4double > f_yi;
      G4int noBinOfPDF;
      G4SimplexDownhill<G4ConvergenceTester>* minimizer;

      G4int noPass;
      G4int noTotal; // Total number of tests

};
#endif

