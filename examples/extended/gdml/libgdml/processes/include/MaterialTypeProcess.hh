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
// $Id: MaterialTypeProcess.hh,v 1.3 2002-06-03 12:09:31 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef MATERIALTYPEPROCESS_H
#define MATERIALTYPEPROCESS_H 1


#include "ProcessingConfigurator.hh"
#include "ProcessingContext.hh"
#include "SAXProcessor.hh"
#include "StateStack.hh"
#include "SAXProcessingState.hh"
#include "SAXStateProcess.hh"
#include "SAXObjectHandle.hh"
#include "SAXComponentFactory.hh"

#include "MaterialType.hh"

#include <cstdlib>
#include <iostream>

class MaterialTypeProcess : public SAXStateProcess, virtual public SAXComponentObject
{
public:
  MaterialTypeProcess( const ProcessingContext* context = 0 )
  : SAXStateProcess( context ), m_obj( 0 )
  {
    // This is not a real process, it's a base class for a real process instead
    // The final inheriting class should set it to a correct value
    //m_tag = "";
  }
  
  virtual ~MaterialTypeProcess()
  {
  }
  
  virtual const SAXComponentObject* Build() const
  {
    return this;
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) { 
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
  }
  
  // Analogical to SAX characters callback, it's called for ignorableWhitespace too!
  virtual void Characters( const std::string& name ) {
  }
  
  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name ) {
    MaterialType* mt_obj = dynamic_cast<MaterialType*>(m_obj);
    
    if( mt_obj == 0 ) {
      // Whoops, fatal problem!
      std::cerr << "MaterialType PROCESS:: bad_cast trouble at: "
                << __FILE__ << ":" << __LINE__ << "\a" << std::endl;
      // Some exception would fit here...
      return;
    }
    
    //std::cout << "MaterialType PROCESS::  : " << name << std::endl;
    
    SAXObject** obj = Context()->GetTopObject();
    
    if( name == "RL" || name == "RLref" ) {      
      mt_obj->set_RLorRLref( name, *obj );
    } else if( name == "AL" || name == "ALref" ) {      
      mt_obj->set_ALorALref( name, *obj );
    } else if( name == "T" || name == "Tref" ) {      
      mt_obj->set_TorTref( name, *obj );
    } else if( name == "P" || name == "Pref" ) {      
      mt_obj->set_PorPref( name, *obj );
    } else {
      // Invalid input
      std::cerr << "MaterialType PROCESS:: invalid tag "
                << name << " on input at: "
                << __FILE__ << ":" << __LINE__ << "\a" << std::endl;
    }
  }
  
  // The name of the state this object will process
//   virtual const std::string& State() const
//   {
//     return m_tag;
//   }

protected:
  // Resulting object, one of isotope, element or material
  SAXObject*  m_obj;
};

// This is going to be a base class, so no factory is needed
//DECLARE_PROCESS_FACTORY(MaterialTypeProcess)

#endif // MATERIALTYPEPROCESS_H
