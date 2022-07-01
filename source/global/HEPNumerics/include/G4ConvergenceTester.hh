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
// G4ConvergenceTester
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
// Positive numbers are assumed for input values

// Author: Tatsumi Koi (SLAC/SCCS)
// --------------------------------------------------------------------
#ifndef G4ConvergenceTester_hh
#define G4ConvergenceTester_hh 1

#include "G4SimplexDownhill.hh"
#include "G4Timer.hh"
#include "globals.hh"

#include <map>
#include <vector>

class G4ConvergenceTester
{
  public:

    G4ConvergenceTester(const G4String& theName = "NONAME");
    ~G4ConvergenceTester();
    G4ConvergenceTester(G4double);

    void AddScore(G4double);

    inline G4ConvergenceTester& operator+=(G4double val)
    {
      this->AddScore(val);
      return *this;
    }

    void ShowHistory(std::ostream& out = G4cout);
    void ShowResult(std::ostream& out = G4cout);
      // Default to G4cout but can be redirected to another ostream

    inline G4double GetValueOfMinimizingFunction(std::vector<G4double> x)
    {
      return slope_fitting_function(x);
    }

    void ComputeStatistics() { calStat(); }
      // Explicitly calculate statistics

    // All accessors check to make sure value is current before returning

    inline G4double GetMean() { CheckIsUpdated(); return mean; }
    inline G4double GetStandardDeviation() { CheckIsUpdated(); return sd; }
    inline G4double GetVariance() { CheckIsUpdated(); return var; }
    inline G4double GetR() { CheckIsUpdated(); return r; }
    inline G4double GetEfficiency() { CheckIsUpdated(); return efficiency; }
    inline G4double GetR2eff() { CheckIsUpdated(); return r2eff; }
    inline G4double GetR2int() { CheckIsUpdated(); return r2int; }
    inline G4double GetShift() { CheckIsUpdated(); return shift; }
    inline G4double GetVOV() { CheckIsUpdated(); return vov; }
    inline G4double GetFOM() { CheckIsUpdated(); return fom; }

  private:

    void calStat();
      // Boolean value of 'statsAreUpdated' is set to TRUE at end of calStat
      // and set to FALSE at end of AddScore().
      // NOTE: A thread lock needs to be put in AddScore() so calStat()
      // is not executed in one thread while AddScore() is adding data

    inline void CheckIsUpdated()
    {
      if(!statsAreUpdated) { calStat(); }
    }

    void calc_grid_point_of_history();
    void calc_stat_history();
    void check_stat_history(std::ostream& out = G4cout);
    G4double calc_Pearson_r(G4int,std::vector<G4double>,std::vector<G4double>);
    G4bool is_monotonically_decrease(const std::vector<G4double>&);
    void calc_slope_fit(const std::vector<G4double>&);
    G4double slope_fitting_function(std::vector<G4double>);

  private:

    G4String name;
    std::map<G4int, G4double> nonzero_histories;  // (ith-history, score value)
    G4int n = 0;        // number of history
    G4double sum = 0.0; // sum of scores;

    G4Timer* timer = nullptr;
    std::vector<G4double> cpu_time;

    G4double mean       = 0.0;
    G4double var        = 0.0;
    G4double sd         = 0.0;
    G4double r          = 0.0;  // relative err sd/mean/sqrt(n)
    G4double efficiency = 0.0;  // rate of non zero score
    G4double r2eff      = 0.0;
    G4double r2int      = 0.0;
    G4double shift      = 0.0;
    G4double vov        = 0.0;
    G4double fom        = 0.0;

    G4double largest             = 0.0;
    G4int largest_score_happened = 0;

    G4double mean_1  = 0.0;
    G4double var_1   = 0.0;
    G4double sd_1    = 0.0;
    G4double r_1     = 0.0;  // relative err sd/mean/sqrt(n)
    G4double shift_1 = 0.0;
    G4double vov_1   = 0.0;
    G4double fom_1   = 0.0;

    G4int noBinOfHistory = 16;
    std::vector<G4int> history_grid;
    std::vector<G4double> mean_history;
    std::vector<G4double> var_history;
    std::vector<G4double> sd_history;
    std::vector<G4double> r_history;
    std::vector<G4double> vov_history;
    std::vector<G4double> fom_history;
    std::vector<G4double> shift_history;
    std::vector<G4double> e_history;
    std::vector<G4double> r2eff_history;
    std::vector<G4double> r2int_history;

    G4double slope = 0.0;
    std::vector<G4double> largest_scores;
    std::vector<G4double> f_xi;
    std::vector<G4double> f_yi;
    G4int noBinOfPDF = 10;
    G4SimplexDownhill<G4ConvergenceTester>* minimizer = nullptr;

    G4int noPass  = 0;
    G4int noTotal = 8;  // Total number of tests

    G4bool statsAreUpdated = true;
    G4bool showHistory     = true;
    G4bool calcSLOPE       = true;
};

#endif
