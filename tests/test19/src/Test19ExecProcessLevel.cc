
#include "Test19ExecProcessLevel.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ios.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleChange.hh"

#include "TstDiscreteProcessReader.hh"
#include "FTFPWrapper.hh"
#include "QGSPWrapper.hh"

void Test19ExecProcessLevel::InitProcess( const TstReader* pset )
{

   // G4String name = (dynamic_cast<const TstDiscreteProcessReader*>(pset))->GetProcessName();
   G4String name = (pset)->GetPhysics();

   ProcessWrapper* pw = 0;
   
   if ( name.find("ftf") != std::string::npos )
   {
      // fProcWrapper = new FTFPWrapper();
      pw = new FTFPWrapper();
   }
   else if ( name.find("qgs") != std::string::npos )
   {
       // fProcWrapper = new QGSPWrapper();
       pw = new QGSPWrapper();
   }
    
   // if (!fProcWrapper) 
   if (!pw) 
   { 
      G4cout 
	     << " generator " << name << " is unavailable"
	     << G4endl;
      exit(1);
   } 
    
   // std::cout << " process = " << fProcWrapper->GetProcessName() << std::endl;
    
   if ( name.find("lund-str-fragm") != std::string::npos )
   {
      // fProcWrapper->UseG4LundStringFragm(true);
      pw->UseG4LundStringFragm(true);
   }
    
   // fProcWrapper->Compose();
   pw->Compose();
   
   fProcWrapper = pw;
   
   return;

}


