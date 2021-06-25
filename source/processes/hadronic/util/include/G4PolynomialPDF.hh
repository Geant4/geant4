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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4PolynomialPDF
//
//      Author:        Jason Detwiler (jasondet@gmail.com)
// 
//      Creation date: Aug 2012
//
//      Description:   Evaluates, generates random numbers from, and evaluates
//      the inverse of a polynomial PDF, its CDF, and its first and second
//      derivative.
//
// -------------------------------------------------------------------

#ifndef G4POLYNOMIALPDF_HH
#define G4POLYNOMIALPDF_HH

#include "globals.hh"
#include <vector>

class G4PolynomialPDF
{
  public:
    G4PolynomialPDF(size_t n = 0, const double* coeffs = nullptr, 
		    G4double x1=0, G4double x2=1);

    ~G4PolynomialPDF();
    // Setters and Getters for coefficients
    inline void SetNCoefficients(size_t n) { fCoefficients.resize(n); fChanged = true; }
    inline size_t GetNCoefficients() const { return fCoefficients.size(); }
    inline void SetCoefficients(const std::vector<G4double>& v) { 
      fCoefficients = v; fChanged = true; Simplify(); 
    }
    inline G4double GetCoefficient(size_t i) const { return fCoefficients[i]; }
    void SetCoefficient(size_t i, G4double value, bool doSimplify);
    void SetCoefficients(size_t n, const G4double* coeffs);
    void Simplify();

    // Set the domain over which random numbers are generated and over which
    // the CDF is evaluated
    void SetDomain(G4double x1, G4double x2);

    // Normalize PDF to 1 over domain fX1 to fX2. Used internally by
    // GetRandomX(), but the user may want to call this as well for evaluation
    // purposes.
    void Normalize();

    // Evaluate (d/dx)^ddxPower f(x) (-1 <= ddxPower <= 2)
    // ddxPower = -1 -> CDF; 
    // ddxPower = 0 -> PDF
    // ddxPower = 1 -> PDF'
    // ddxPower = 2 -> PDF''
    G4double Evaluate(G4double x, G4int ddxPower = 0);

    // Generate a random number from this PDF
    G4double GetRandomX();

    // Set the tolerance to within negative minima are checked
    inline void SetTolerance(G4double tolerance) { fTolerance = tolerance; }

    // Find a value x between x1 and x2 at which ddxPower[PDF](x) = p.
    // ddxPower = -1 -> CDF; 
    // ddxPower = 0 -> PDF
    // ddxPower = 1 -> PDF'
    // (ddxPower = 2 not implemented)
    // Solves analytically when possible, and otherwise uses the Newton-Raphson
    // method to find the zero of ddxPower[PDF](x) - p.
    // If not found in range, returns the nearest boundary.
    // Beware that if x1 and x2 are not set carefully there may be multiple
    // solutions, and care is not taken to select a particular one among them.
    // Returns x2 on error
    G4double GetX( G4double p, G4double x1, G4double x2, G4int ddxPower = 0, 
                   G4double guess = 1.e99, G4bool bisect = true );
    inline G4double EvalInverseCDF(G4double p) { return GetX(p, fX1, fX2, -1, fX1 + p*(fX2-fX1)); }
    G4double Bisect( G4double p, G4double x1, G4double x2 );

    void Dump();

  protected:
    // Checks for negative values between x1 and x2. Used by GetRandomX()
    G4bool HasNegativeMinimum(G4double x1, G4double x2);

    G4double fX1;
    G4double fX2;
    std::vector<G4double> fCoefficients;
    G4bool   fChanged;
    G4double fTolerance;
    G4int    fVerbose;
};

#endif
