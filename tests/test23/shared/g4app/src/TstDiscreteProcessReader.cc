
#include "TstDiscreteProcessReader.hh"

void TstDiscreteProcessReader::ProcessLine( G4String line )
{
   
   TstReader::ProcessLine( line );
   
   if ( line != "#generator" ) return;
   
   if(line == "#generator") 
   {
        (*fInStream) >> fPhysics;
   }

   return;

}
