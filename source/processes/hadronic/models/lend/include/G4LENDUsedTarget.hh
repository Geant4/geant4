#ifndef G4LENDUsedTarget_h
#define G4LENDUsedTarget_h 1

// Class Description
// Container of LEND (Low Energy Nuclear Data) target (nucleus) 
// LEND is Geant4 interface for GIDI (General Interaction Data Interface) 
// which gives a discription of nuclear and atomic reactions, such as
//    Binary collision cross sections
//    Particle number multiplicity distributions of reaction products
//    Energy and angular distributions of reaction products
//    Derived calculational constants
// GIDI is developped at Lawrence Livermore National Laboratory
// Class Description - End

// 071025 First implementation done by T. Koi (SLAC/SCCS)
// 101118 Name modifications for release T. Koi (SLAC/PPA)

#include "G4LENDHeader.hh"
#include "G4ParticleDefinition.hh"

class G4LENDUsedTarget 
{

   public:

      G4LENDUsedTarget( G4ParticleDefinition* pd , G4String Evaluation , G4int iZ , G4int iA , G4int iM = 0 )
      : allow_nat ( false ) 
      , allow_anything ( false ) 
      , min_Z ( 0 )
      , max_Z ( 113 ) 
      , min_A ( 1 )
      , max_A ( 278 )
      , min_M ( 0 )
      , max_M ( 10 )
      {

         proj = pd;

         wanted_Z = iZ;
         wanted_A = iA;
         wanted_M = iM;
         wanted_Evaluation = Evaluation;

         actual_Z = -1;
         actual_A = -1;
         actual_M = -1;
         actual_Evaluation = "na";

         searchTarget();
      }

      ~G4LENDUsedTarget(){;};

      void AllowNat()
      {
         allow_nat = true;
         searchTarget();
      }; 

      void AllowAny()
      {
         allow_anything = true;
         searchTarget();
      }; 

      G4int GetWantedZ(){ return wanted_Z; };
      G4int GetWantedA(){ return wanted_A; };
      G4int GetWantedM(){ return wanted_M; };

      G4int GetActualZ(){ return actual_Z; };
      G4int GetActualA(){ return actual_A; };
      G4int GetActualM(){ return actual_M; };

      G4String GetWantedEvaluation(){ return wanted_Evaluation; };
      G4String GetActualEvaluation(){ return actual_Evaluation; };

      G4GIDI_target* GetTarget(){ return target; };

   private:

      void searchTarget();

      G4ParticleDefinition* proj;

      G4int wanted_Z;
      G4int wanted_A;
      G4int wanted_M;

      G4String wanted_Evaluation;

      G4bool allow_nat; 
      G4bool allow_anything; 

      G4GIDI_target* target;
      
      G4int actual_Z;
      G4int actual_A;
      G4int actual_M;
      G4String actual_Evaluation;

      G4int min_Z;
      G4int max_Z;
      G4int min_A;
      G4int max_A;
      G4int min_M;
      G4int max_M;
};

#endif
