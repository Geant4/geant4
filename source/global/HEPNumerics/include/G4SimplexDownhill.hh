//
// Minimization of a functionof n variables
//
// Reference 
// "A Simplex method for function minimization" By J. A. Nelder and R. Mead, Computer Journal, 7, 308 (1965)
// and also 
// "Numerical Recipes in C: the art of scientific computing" By William H., Cambridge University Press ISBN 0521437202 (1992)
//
// Koi, Tatsumi (SLAC/SCCS)
// 

#ifndef G4SimplexDownhill
#define G4SimplexDownhill_h

#include "globals.hh"

#include <vector> 
#include <algorithm> 

template<class T> class G4SimplexDownhill
{

      T* target;
   public: 

      G4SimplexDownhill( T* tp , G4int n ){ target = tp; init(); numberOfVariable = n; };

      ~G4SimplexDownhill();

      G4double GetMinimum();

      std::vector< G4double > GetMinimumPoint();


   private:
      G4double getValue( std::vector< G4double > x ) { return target->GetValueOfMinimizingFunction( x ); };

      G4int numberOfVariable; 

      void initialize();
      std::vector< std::vector< G4double > > currentSimplex;

      void calHeights();
      std::vector< G4double > currentHeights;
      G4double currentValue;

      std::vector< G4double > calCentroid( G4int );

      G4bool isItGoodEnough();

      std::vector< G4double > getReflectionPoint( std::vector< G4double > , std::vector< G4double > );
      std::vector< G4double > getExpansionPoint( std::vector< G4double > , std::vector< G4double > );
      std::vector< G4double > getContractionPoint( std::vector< G4double > , std::vector< G4double > );

      void doDownhill();

      void init();

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
