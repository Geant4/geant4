
#include "TstPhysListReader.hh"

void TstPhysListReader::ProcessLine( G4String line )
{
   
   TstReader::ProcessLine( line );
   
   if ( line != "#physicslist" ) return;
   
   if(line == "#physicslist") 
   {
        (*fInStream) >> fPhysics;
   }

   return;

}
