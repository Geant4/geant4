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
//      File name:     G4LegendrePolynomial
//
//      Author:        Jason Detwiler (jasondet@gmail.com)
// 
//      Creation date: February 2015
//
//      Modifications: 
//      
//      Legendre Polynomial
//
// -------------------------------------------------------------------

#ifndef G4LEGENDREPOLYNOMIAL_HH
#define G4LEGENDREPOLYNOMIAL_HH

#include "globals.hh"
#include <vector>
#include <map>

class G4LegendrePolynomial
{
  public:
    // Access to coefficients
    static size_t GetNCoefficients(size_t order) { return order+1; }
    G4double GetCoefficient(size_t i, size_t order);

    // Evaluation functions
    G4double EvalLegendrePoly(G4int order, G4double x);
    G4double EvalAssocLegendrePoly(G4int l, G4int m, G4double x,
                                   std::map<G4int, std::map<G4int, G4double> >* cache = NULL);

  protected: // Cache coefficients for speed
    void BuildUpToOrder(size_t order);
    std::vector< std::vector<G4double> > fCoefficients;
};

#endif
