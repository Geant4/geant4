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
// $Id: paraSubscriber.cc,v 1.3 2002-06-03 12:09:35 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#include "SAXSubscriber.hh"
#include "SAXComponentFactory.hh"

#include "GDMLProcessor.hh"
#include "GDMLExpressionEvaluator.hh"

#include "para.hh"

#include "G4VSolid.hh"
#include "G4Para.hh"

#include <iostream>

class paraSubscriber : public SAXSubscriber
{
public:
  virtual const SAXComponentObject* Build() const {
    return this;
  }

public:
  paraSubscriber() {
    Subscribe( "para" );
  }
  virtual ~paraSubscriber() {
  }
   
  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object ) {
    //std::cout << "PARA SUBSCRIBER:: ";
    if( object != 0 ) {
      try {
        const para* obj = dynamic_cast<const para*>( object );    
        
        //std::cout << "GOT PARA " << obj->get_name() << std::endl;
      
        GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
      
        const std::string& name = obj->get_name();
        std::string lunit = obj->get_lunit();
        std::string aunit = obj->get_aunit();
        
        std::string sval = obj->get_x();
        sval += "*"+lunit;
        double dx = calc->Eval( sval ); dx = dx/2.;
        sval = obj->get_y();
        sval += "*"+lunit;
        double dy = calc->Eval( sval ); dy = dy/2.;
        sval = obj->get_z();
        sval += "*"+lunit;
        double dz = calc->Eval( sval ); dz = dz/2.;
        sval = obj->get_alpha();
        sval += "*"+aunit;
        double dalpha = calc->Eval( sval );
        sval = obj->get_theta();
        sval += "*"+aunit;
        double dtheta = calc->Eval( sval );
        sval = obj->get_phi();
        sval += "*"+aunit;
        double dphi = calc->Eval( sval );
        
//         std::cout << "x:     " << obj->get_x()     << lunit << " dx:     " << dx     << std::endl;
//         std::cout << "y:     " << obj->get_y()     << lunit << " dy:     " << dy     << std::endl;
//         std::cout << "z:     " << obj->get_z()     << lunit << " dx:     " << dz     << std::endl;
//         std::cout << "alpha: " << obj->get_alpha() << aunit << " dalpha: " << dalpha << std::endl;
//         std::cout << "theta: " << obj->get_theta() << aunit << " dtheta: " << dtheta << std::endl;
//         std::cout << "phi:   " << obj->get_phi()   << aunit << " dphi:   " << dphi   << std::endl;

        G4VSolid* newobj = new G4Para( name, dx, dy, dz, dalpha, dtheta, dphi);
      
        GDMLProcessor::GetInstance()->AddSolid( name, newobj );      
      } catch(...) {
        std::cerr << "GOT INTO BAD_CAST TROUBLE!" << std::endl;
      }
    } else {
      std::cerr << "GOT ZERO DATA POINTER!" << std::endl;
    }
    delete object;
  }
};

DECLARE_SUBSCRIBER_FACTORY(paraSubscriber)

