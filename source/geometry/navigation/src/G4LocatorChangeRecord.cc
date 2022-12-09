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
// G4LocatorChangeRecord class implementation
//
// Author: John Apostolakis, 28.08.19 - First version
// --------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <locale>
// #include <cassert>

#include "G4LocatorChangeRecord.hh"

// Must correspond exactly with the count of the enum EChangeLocation
//
const char * G4LocatorChangeRecord::fNameChangeLocation[] =
  { "Invalid", "Unknown", "Initialising", "IntersectsAF", "IntersectsFB",  // 5
    "NoIntersections-AForFB",  "RecalculatedB",         // 2
    "InsertingMidPoint", "RecalculatedB-2ndHalf",       // 2
    "Level Pop" };

// --------------------------------------------------------------------
// 
std::ostream& G4LocatorChangeRecord::ReportVector ( std::ostream& os,
                                              const std::string & name,
                       const std::vector<G4LocatorChangeRecord> & vecRec )
{
  using std::setw;
  G4int prec= 16;
  if( vecRec.size() == 0 )
  {
     os << "Locator Change Record for " << name << "  is empty" << G4endl;
     return os;
  }
  
  G4long oldprc = os.precision(prec);

  // std::vector<G4LocatorChangeRecord>::const_iterator
  auto itRec
     = std::vector<G4LocatorChangeRecord>::const_iterator(vecRec.cbegin());

  os << setw( 7 ) << "Change#"  << "  " 
     << setw( 4 ) << "Iter"  << "  "
     << std::left
     << setw( prec+9 ) << "Length" << "  "
     << setw( 15 ) << "Code-Location"  << "  "
     << G4endl;
  os << "====================================================================="
     << G4endl;
  
  do
  {
     auto locationCode= (*itRec).GetLocation();
     os << std::internal
        << setw( 7 ) << (*itRec).GetCount() << "  "       // Event Count
        << setw( 4 ) << (*itRec).GetIteration()  << "  "
        << std::left
        << setw( prec+9 ) << (*itRec).GetLength() << "  "
        << setw( 2  ) << locationCode  << "  "  // location enum
        << setw( 15 ) << fNameChangeLocation[ locationCode ]
        << std::internal
        ;
     os << G4endl;
     ++itRec;

  } while ( itRec != vecRec.cend() );

  os.precision(oldprc);
  return os;
}

// --------------------------------------------------------------------
// 
std::ostream& 
G4LocatorChangeRecord::ReportEndChanges (
   std::ostream& os,
   const std::vector<G4LocatorChangeRecord> & startA,
   const std::vector<G4LocatorChangeRecord> &   endB )
{
  using std::setw;   
  G4int prec= 16;
  const G4bool confirm = true;
  G4long oldprc = os.precision(prec);

  std::vector<G4LocatorChangeRecord>::const_iterator itrecA, itrecB;
  itrecA= startA.begin();
  itrecB= endB.begin();
  os << G4endl;
  os << "=========================================================================================";
  os << G4endl << " ** Change records: " << G4endl;
  os << " *     endPoints  A (start) and B (end): combined changes of AB intervals" << G4endl;
  os << " *     Sizes of change records:  start(A) : " << startA.size() 
     << "  end(B) : " <<   endB.size() << G4endl;
  os << "========================================================================================="    
     << G4endl;

  os << setw( 7 ) << "Change#"  << "  "
     << setw( 4 ) << "Iter"  << "  "
     << setw( 20 ) << "CodeLocation"  << "  "
     << setw( prec+9 ) << "Length-A (start)" << "  " 
     << setw( prec+9 ) << "Length-B (end)" << "  "
     << G4endl;  
  os << "=====================================================================";

  
  auto eventA = (*itrecA).GetCount();
  auto eventB = (*itrecB).GetCount();

  G4bool isLastA= false;
  G4bool isLastB= false;

  G4int maxEvent = std::max( startA[ startA.size() - 1 ].GetCount() ,
                             endB[   endB.size() - 1 ].GetCount() );
  G4int prevA = -1;
  G4int prevB = -1;  

  G4bool advanceA= false, advanceB= false;
  do
  {
     advanceA= false;
     advanceB= false;

     if( ((G4int)eventA>prevA) && ((G4int)eventB>prevB) )
     {
        auto codeLocA= (*itrecA).GetLocation();
        
        os << G4endl;                   
        os << setw( 7 ) << eventA  << "  " 
           << setw( 4 ) << (*itrecA).GetIteration()  << "  "
           << setw( 3 ) << codeLocA << " " 
           << setw( 15 ) << fNameChangeLocation[ codeLocA ] << " " 
           << setw( prec+9 ) << (*itrecA).GetLength() << "  " 
           << setw( prec+9 ) << (*itrecB).GetLength() << "  ";
        if( confirm )
        { 
          os << setw( 4 ) << (*itrecB).GetIteration()  << "  "
             << setw( 15 ) << (*itrecB).GetLocation();
        }
     }
     else
     {
        if ( (G4int)eventA > prevA )
        {
           auto codeLocA = (*itrecA).GetLocation();           
           os << G4endl;           
           os << setw( 7 ) << (*itrecA).GetCount() << "  " 
              << setw( 4 ) << (*itrecA).GetIteration()  << "  "
              << setw( 3 ) << codeLocA << " " 
              << setw( 15 ) << fNameChangeLocation[ codeLocA ] << " "            
              << setw( prec+9 ) << (*itrecA).GetLength() << "  " 
              << setw( prec+9 ) << "       " << "  ";
        }
        else
        {
           // assert( (G4int)eventB > prevB );
           auto codeLocB = (*itrecB).GetLocation();
           
           os << G4endl;           
           os << setw( 7 ) << eventB  << "  " 
              << setw( 4 ) << (*itrecB).GetIteration()  << "  "
              << setw( 3 ) << codeLocB << " " 
              << setw( 15 ) << fNameChangeLocation[ codeLocB ] << " "                          
              << setw( prec+9 ) << "       " << "  "
              << setw( prec+9 ) << (*itrecB).GetLength() << "  " ;           
        }
     }

     prevA= eventA;
     prevB= eventB;
     
     auto nextA= itrecA;
     auto nextB= itrecB;
     
     G4int nextAct = maxEvent, nextBct = maxEvent;
     ++nextA;
     ++nextB;
     if ( nextA != startA.end() ) { nextAct = (*nextA).GetCount(); }
     if ( nextB !=   endB.end() ) { nextBct = (*nextB).GetCount(); }

     isLastA=  ( nextA >= startA.end() );
     isLastB=  ( nextB >=   endB.end() );
     
     advanceA= ( nextAct <= nextBct ) && !isLastA;
     advanceB= ( nextBct <= nextAct ) && !isLastB;

     if( advanceA )
     {
        ++itrecA;
        if( !isLastA ) { eventA = (*itrecA).GetCount(); }
        else { eventA = maxEvent; }
     }
     
     if( advanceB )
     {
        ++itrecB;
        if( !isLastB ) { eventB = (*itrecB).GetCount(); } 
        else { eventB = maxEvent; }
     }

     // Checks
     if( isLastA !=  ( nextA == startA.end() ) )
     {
        os << G4endl;
        os << "  Checking isLastA= " << isLastA << " vs expected :  "
           << ( itrecA == startA.end() );
        os << " BAD --- ERROR " << G4endl;
     }
     if( isLastB !=  ( nextB == endB.end() ) )
     {
        os << G4endl;        
        os << "  Checking isLastB= " << isLastB << " vs expected :  "
           << ( itrecB ==  endB.end() );        
        os << " BAD --- ERROR " << G4endl;        
     }
       
  } while ( ! ( isLastA && isLastB ) );
  
  os << G4endl;           
  os.precision(oldprc);
  return os;
}

// --------------------------------------------------------------------
// Streaming operator dumping record 
//
std::ostream& operator<< ( std::ostream& os, const G4LocatorChangeRecord& e )
{
  return e.StreamInfo(os);
}
      
// --------------------------------------------------------------------
// Stream object contents to an output stream
//
std::ostream& G4LocatorChangeRecord::StreamInfo(std::ostream& os) const
{
  G4long oldprc = os.precision(16);
  os << "  count = " << fEventCount
     << "  iter= "   << fIteration
     << "  Location code = " << fCodeLocation
     << "  Length = " << GetLength() << G4endl;
  os.precision(oldprc);
  return os;
}

// --------------------------------------------------------------------
// Streaming operator
//
std::ostream& operator << ( std::ostream& os,
                            const std::vector<G4LocatorChangeRecord> & vecR )
{
  G4LocatorChangeRecord::ReportVector( os, "", vecR );
  return os;
}

// --------------------------------------------------------------------
// 
const char *G4LocatorChangeRecord::GetNameChangeLocation( EChangeLocation loc )
{
  return fNameChangeLocation[loc];
}
