//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: BooleanSolidTypeProcess.hh,v 1.3 2002-06-03 12:09:31 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef BOOLEANSOLIDTYPEPROCESS_H
#define BOOLEANSOLIDTYPEPROCESS_H 1

#include "BooleanSolidType.hh"

#include "SolidTypeProcess.hh"

class BooleanSolidTypeProcess : public SolidTypeProcess
{
public:
  BooleanSolidTypeProcess( const ProcessingContext* context = 0 )
  : SolidTypeProcess( context ) {
  }
  
  virtual ~BooleanSolidTypeProcess() {
  }
  
  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name ) {
    
    SAXObject** so = Context()->GetTopObject();
    BooleanSolidType* bsobj = dynamic_cast<BooleanSolidType*>( m_obj );
    
    if( name == "first" || name == "second" ) {
      // Add the first and second elements into content model
      bsobj->add_content( name, *so );
    } else {
      // Add either choice of position or rotation elements and their references
      ContentChoice* choice        = new ContentChoice;
      ContentGroup::ContentItem ci = { name, *so };
      choice->set_content( ci );
      //std::cout << "BOOLEAN SOLID PROCESS adding choice " << name << std::endl;
      bsobj->add_content( "choice", choice );
    }
    
    SolidTypeProcess::StackPopNotify( name );
  }
};

#endif // BOOLEANSOLIDTYPEPROCESS_H
