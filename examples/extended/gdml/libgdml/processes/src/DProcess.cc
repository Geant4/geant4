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
// $Id: DProcess.cc,v 1.2 2002-06-03 12:09:32 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#include "QuantityTypeProcess.hh"

#include "DorDref.hh"

class DProcess : public QuantityTypeProcess
{
public:
  DProcess( const ProcessingContext* context = 0 )
  : QuantityTypeProcess( context ) {
  }
  
  virtual ~DProcess() {
  }
  
  // Analogical to SAX startElement callback
  virtual void StartElement( const std::string& name, const ASCIIAttributeList& attrs ) {    

    SAXObject** obj = Context()->GetTopObject();
    *obj = new D;
    m_obj = *obj;

    QuantityTypeProcess::StartElement( name, attrs );    
  }
  
  // Analogical to SAX endElement callback
  virtual void EndElement( const std::string& name )
  {
    QuantityTypeProcess::EndElement( name );
  }
  
  // The name of the state this object will process
  virtual const std::string& State() const
  {
    static std::string m_tag = "D";
    return m_tag;
  }
};

DECLARE_PROCESS_FACTORY(DProcess)

