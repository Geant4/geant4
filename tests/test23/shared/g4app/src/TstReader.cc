
#include "TstReader.hh"

#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"

TstReader::TstReader()
   : fInStream(0), fEndConfig(false),
     fNEvents(1000), 
     fBeamPart("proton"),
     fBeamEnergy(0.),
     fBeamMomentum(1000.*MeV),
     fDirection( G4ThreeVector(0.,0.,1.) ),
     fPosition( G4ThreeVector(0.,0.,0.) ),
     fTime(0.),
     fTargetMaterial("C"),
     fTargetSize( G4ThreeVector(100., 100., 100. ) ),
     fTargetShape("G4Box"),
     fPhysics(""),
     fStep(0.01*micrometer),
     fRndmSeed(135799753),
     fJobID(-1),
     fVerbose(-1),
     fExpDataSet("None"),
     fForseResDecay(false)    
{
   
   fInStream = new std::ifstream();

}

TstReader::~TstReader()
{

   if ( !fInStream )
   {
      if ( fInStream->is_open() ) fInStream->close();
      delete fInStream;
   }

}

void TstReader::Help() 
{

   G4cout << "Available configuration commands are: " << G4endl;
   G4cout << "#events" << G4endl;
   G4cout << "#particle" << G4endl;
   G4cout << "#energy(MeV)" << G4endl;
   G4cout << "#momentum(MeV/c)" << G4endl;
   G4cout << "#step(mm)" << G4endl;
   G4cout << "#material" << G4endl;
   G4cout << "#target-geom(mm)" << G4endl;
   G4cout << "#generator" << G4endl;
   G4cout << "#verbose" << G4endl;
   G4cout << "#position(mm)" << G4endl;
   G4cout << "#direction" << G4endl;
   G4cout << "#time(ns)" << G4endl; // why would I need this ?
   G4cout << "isHARP" << G4endl;
   G4cout << "isNA61" << G4endl;
   G4cout << "isNA49" << G4endl;
   // G4cout << "isMIPP" << G4endl;
   // G4cout << "isITEP" << G4endl;
   // G4cout << "isBNL" << G4endl;
   G4cout << "#forseResonanceDecay" << G4endl;
  
//
// for parallel processing
//
   G4cout << "#randomSeed" << G4endl;
   G4cout << "#jobID" << G4endl;

   G4cout << "#run" << G4endl;
   G4cout << "#exit" << G4endl;

   return;

}

void TstReader::OpenAppConfig( std::string conf )
{

   fInStream->open( conf.c_str() );
   if ( !fInStream->is_open() )
   {
    G4cout << "Input file <" << conf << "> does not exist! Exit" << G4endl;
    exit(1);
   }   

   return;

}

void TstReader::CloseAppConfig()
{
   
   if ( fInStream->is_open() ) fInStream->close();
   
   return;
   
}

void TstReader::ProcessConfig()
{

   G4String line;
   
   do
   {
      (*fInStream) >> line;
      if ( line == "#exit" )
      {
         fEndConfig = true;
	 break;
      }
      else if ( line == "#run" )
      {
         break;
      }
      else
      {
         ProcessLine( line );
      }
   
   } while(!fEndConfig);

   return;

}

void TstReader::ProcessLine( G4String line )
{

      std::string line1;
      G4double nx, ny, nz;
      
      if(line == "#particle") 
      {
        (*fInStream) >> fBeamPart;
      } 
      else if(line == "#momentum(MeV/c)") 
      {
        (*fInStream) >> fBeamMomentum;
        fBeamMomentum *= MeV;
      } 
      else if(line == "#events") 
      {
        (*fInStream) >> line1;
        std::istringstream is(line1);
        is >> fNEvents;
        // G4cout << "nevt : " << nevt << G4endl;
      } 
      else if(line == "#step(mm)") 
      {
        (*fInStream) >> fStep;
        fStep *= mm;
      } 
      else if(line == "#material") 
      {
        (*fInStream) >> fTargetMaterial;
      } 
      else if ( line == "#target-geom(mm)" )
      {
         (*fInStream) >> nx >> ny >> nz >> fTargetShape;
	 fTargetSize = G4ThreeVector(nx*mm, ny*mm, nz*mm);
      }
      else if(line == "#verbose") 
      {
        (*fInStream) >> fVerbose;
      } 
      else if(line == "#position(mm)") 
      {
	(*fInStream) >> nx >> ny >> nz;
        fPosition = G4ThreeVector(nx*mm, ny*mm, nz*mm);
      } 
      else if(line == "#direction") 
      {
        (*fInStream) >> nx >> ny >> nz;
        if( (nx*nx+ny*ny+nz*nz) > 0.0) 
	{
          fDirection = G4ThreeVector(nx, ny, nz);
          fDirection = fDirection.unit();
	}
      } 
      else if(line == "#time(ns)") 
      {
        (*fInStream) >> fTime;
        fTime *= ns;
      } 
//
// Exp. Data Set type
//
      else if ( line == "#isNA49" )
      {
         SetExpDataSet( "NA49" );
      }
      else if ( line == "#isNA61" )
      {
         SetExpDataSet( "NA61" );
      } 
      else if ( line == "#isHARP" )
      {
         SetExpDataSet( "HARP" );
      }
      else if ( line == "#isITEP" )
      {
         SetExpDataSet( "ITEP" );
      }
      else if ( line == "#isBNL" )
      {
         SetExpDataSet( "BNL" );
      }
      else if ( line == "#isMIPP" )
      {
         SetExpDataSet( "MIPP" );
      }
//
// needed for parallel processing
//
      else if ( line == "#randomSeed" )
      {
        (*fInStream) >> fRndmSeed;
      }
      else if ( line == "#jobID" )
      {
	(*fInStream) >> fJobID ;
      }
      else if ( line == "#forseResonanceDecay" )
      {
         fForseResDecay = true;
      }

   return;

}
