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
// $Id: volumeProcess.cc,v 1.3 2002-06-03 12:09:32 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef SINGLEPLACEMENTTYPEPROCESS_H
#define SINGLEPLACEMENTTYPEPROCESS_H 1

#include "ProcessingConfigurator.hh"
#include "ProcessingContext.hh"
#include "SAXProcessor.hh"
#include "StateStack.hh"
#include "SAXProcessingState.hh"
#include "SAXStateProcess.hh"
#include "SAXComponentFactory.hh"

#include "volume.hh"

class volumeProcess : public SAXStateProcess, virtual public SAXComponentObject
{
public:
  volumeProcess( const ProcessingContext* context = 0 )
  : SAXStateProcess( context ), m_obj( 0 ) {
  }
  
  virtual ~volumeProcess() {
  }
  
  virtual const SAXComponentObject* Build() const {
    return this;
  }

  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {
    //std::cout << "PROCESS::START OF TAG  : " << name << std::endl;
    SAXObject** obj = Context()->GetTopObject();
    
    volume* vo = new volume;
    
    vo->set_name( attrs.getValue( "name" ) );
    
    m_obj = vo;
    *obj  = vo;
  }

  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name ) {
    //std::cout << "PROCESS::END OF TAG  : " << name << " " << std::endl;
  }

  // Analogical to SAX characters callback, it's called for ignorableWhitespace too!
  virtual void Characters( const std::string& name ) {
  }
  
  // Invoked whenever one of the daughter state processes has been popped-out of the state stack
  // The name passed-in as the argument is the name of the XML element for which that's been done
  virtual void StackPopNotify( const std::string& name ) {
    //std::cout << "PROCESS::volume NOTIFIED AFTER THE TAG: " << name << std::endl;
    SAXObject** so = Context()->GetTopObject();
    volume* vobj = dynamic_cast<volume*>( m_obj );
    vobj->add_content( name, *so );
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const
  {
    static std::string tag = "volume";
    return tag;
  }

private:
  SAXObject* m_obj;
};

DECLARE_PROCESS_FACTORY(volumeProcess)
    
#endif // SINGLEPLACEMENTTYPEPROCESS_H
