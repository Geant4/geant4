// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:    G4QMDMeanField.hh 
//
//      Author: Koi, Tatsumi (tkoi@slac.stanford.edu)       
// 
//      Creation date: 29 March 2007
// -----------------------------------------------------------------------------

#ifndef G4QMDMeanField_hh
#define G4QMDMeanField_hh

#include "G4QMDSystem.hh"
#include "G4QMDNucleus.hh"

class G4QMDMeanField 
{
   public:
      G4QMDMeanField();
      ~G4QMDMeanField();

      void SetSystem ( G4QMDSystem* aSystem );
      void SetNucleus ( G4QMDNucleus* aSystem );
      G4QMDSystem* GetSystem () {return system; };

      void Cal2BodyQuantities(); 
      void Cal2BodyQuantities( G4int ); 

      void CalGraduate(); 

      G4bool IsPauliBlocked( G4int ); 

      G4double GetTotalPotential(); 

      void DoPropagation( G4double );

      std::vector< G4QMDNucleus* > DoClusterJudgment();

      G4double GetRR2( G4int i , G4int j ) { return rr2[i][j]; };

      G4double GetRHA( G4int i , G4int j ) { return rha[i][j]; };
      G4double GetRHE( G4int i , G4int j ) { return rhe[i][j]; };
      G4ThreeVector GetFFr( G4int i ) { return ffr[i]; };
      G4ThreeVector GetFFp( G4int i ) { return ffp[i]; };

      std::vector< G4double > GetLocalDensity();  
      std::vector< G4double > GetDepthOfPotential();  

   private:
      G4double calPauliBlockingFactor( G4int ); 

      G4QMDSystem* system; 

      G4double rclds;

      G4double hbc , rho0;
      G4double epsx , epscl; 

      G4double cpc;

      G4int icoul, irelcr;
      G4double gamm, c0, c3, cs, cl, wl; 
      G4double c0w, c3w, clw, c0sw;

      G4double c0g,c3g,csg,pag; 

      G4double cpw,cph;
       
      // 2 Body Quantities 
      std::vector < std::vector < G4double > > rr2;    
      std::vector < std::vector < G4double > > pp2;    
      std::vector < std::vector < G4double > > rbij;    

      // Gauss 
      std::vector < std::vector < G4double > > rha;    

      // Coulomb
      std::vector < std::vector < G4double > > rhe;    
      std::vector < std::vector < G4double > > rhc;    
                                         
      std::vector < G4ThreeVector > ffr;    
      std::vector < G4ThreeVector > ffp;    
      std::vector < G4double > rh3d;    
       
};

#endif
