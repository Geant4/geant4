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
//
// ----------------------------------------------------------------------
// Class G4StatAnalysis
//
// Class description:
//
//  Class for statistical analysis of random variable
//
//  Adapted
//  Lux, I.
//      Monte Carlo particle transport methods: neutron and photon
//      calculations/authors, Ivan Lux and Laszlo Koblinger.
//      ISBN 0-8493-6074-9
//  1. Neutron transport theory. 2. Photon transport theory.
//  3. Monte Carlo method. I. Koblinger, Laszlo. II. Title.
//  QC793.5.N4628L88 1990
//  530.1 '38—dc20
//
// https://gnssn.iaea.org/NSNI/Shared%20Documents/OPEN%20Shared%20Files/MonteCarloParticleTransportMethodsNeutronAndPhotonCalculations.pdf
//
//

#ifndef G4StatAnalysis_hh_
#define G4StatAnalysis_hh_

//----------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>
#include <limits>
#include <fstream>
#include <cmath>

#include "globals.hh"
#include "tls.hh"

#include "G4Types.hh"
#include "G4Timer.hh"
#include "G4ios.hh"
#include "G4Allocator.hh"

class G4StatAnalysis
{
public:
    inline G4StatAnalysis();
    inline ~G4StatAnalysis() { }

public:
    // Accumulated values
    inline G4double         GetMean() const;
    inline const G4double&  GetSum() const;
    inline const G4double&  GetSumSquared() const;
    inline const G4double&  GetSum1() const;
    inline const G4double&  GetSum2() const;
    inline const G4int&     GetHits() const;
    inline G4int            GetNumNonZero() const;
    inline G4int            GetNumZero() const;

    // Some control over accumulated variables
    inline void SetSum(const G4double& val);
    inline void SetSumSquared(const G4double& val);
    inline void SetSum1(const G4double& val);
    inline void SetSum2(const G4double& val);
    inline void SetHits(const G4int& val);
    inline void SetZero(const G4int& val);

    // Computed values
    inline G4double GetFOM() const;
    inline G4double GetRelativeError() const;
    inline G4double GetStdDev() const;
    inline G4double GetVariance() const;
    inline G4double GetCoeffVariation() const;
    inline G4double GetEfficiency() const;
    inline G4double GetR2Int() const;
    inline G4double GetR2Eff() const;

    // Conversion
    inline operator G4double() const;

    // Modifications
    inline void Reset();
    inline void Add(const G4double& _val, const G4double& _weight = 1.0);
    inline void Rescale(const G4double& factor);

    // Output
    inline void PrintInfo(std::ostream& os, const std::string& = "") const;

    // Operators
    inline G4StatAnalysis& operator+=(const G4double& _val);
    inline G4StatAnalysis& operator/=(const G4double& _val);
    inline G4StatAnalysis& operator+=(const G4StatAnalysis&);
    inline G4StatAnalysis& operator-=(const G4StatAnalysis&);

    // Allocators
    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // Timing (member functions)
    inline G4double GetCpuTime() const;
    // Timing (static functions)
    static tms*& GetCpuClock()
    {
        G4ThreadLocalStatic tms* _instance = nullptr;
        if(!_instance)
        {
            _instance = new tms;
            times(_instance);
        }
        return _instance;
    }
    // Note: this above implementation was implemented in such a way as to
    // conserve memory by eliminated every instance from requiring their own
    // timing variables. The ResetCpuClock function below is called at the
    // beginning of the run (G4Run constructor) to attempt to ensure the
    // FOM is not skewed by multiple runs -- it may be necessary to
    // manually invoke in some situations
    static void ResetCpuClock()
    {
        tms*& _clock = GetCpuClock();
        times(_clock);
    }

private:
    // summation of each history^1
    G4double    fSum1;
    // summation from each history^2
    G4double    fSum2;
    // number of scoring histories
    G4int       fHits;
    // number of histories that were not greater than 0.0
    G4int       fZero;

public:
    // friend operator for output
    friend std::ostream& operator<<(std::ostream& os, const G4StatAnalysis& obj)
    {
        obj.PrintInfo(os);
        return os;
    }
    // friend operator for addition
    friend const G4StatAnalysis operator+(const G4StatAnalysis& lhs,
                                          const G4StatAnalysis& rhs)
    {
        return G4StatAnalysis(lhs) += rhs;
    }
    // friend operator for subtraction
    friend const G4StatAnalysis operator-(const G4StatAnalysis& lhs,
                                          const G4StatAnalysis& rhs)
    {
        return G4StatAnalysis(lhs) -= rhs;
    }
};

#include "G4StatAnalysis.icc"

#endif
