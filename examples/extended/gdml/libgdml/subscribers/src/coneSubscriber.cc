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
// $Id: coneSubscriber.cc,v 1.3 2002-06-03 12:09:35 radoone Exp $
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

#include "cone.hh"

#include "G4VSolid.hh"
#include "G4Cons.hh"

#include <iostream>

class coneSubscriber : public SAXSubscriber
{
public:
  virtual const SAXComponentObject* Build() const {
    return this;
  }

public:
  coneSubscriber() {
    Subscribe( "cone" );
  }
  virtual ~coneSubscriber() {
  }
   
  // The activation callback invoked by SAXProcessor whenever it has
  // a new object created from XML and a corresponding subcriber exists
  virtual void Activate( const SAXObject* object ) {
    //std::cout << "CONE SUBSCRIBER:: ";
    if( object != 0 ) {
      try {
        const cone* obj = dynamic_cast<const cone*>( object );    
        
        //std::cout << "GOT CONE " << obj->get_name() << std::endl;
      
        GDMLExpressionEvaluator* calc = GDMLProcessor::GetInstance()->GetEvaluator();
      
        std::string lunit = obj->get_lunit();
        std::string aunit = obj->get_aunit();
        const std::string& name = obj->get_name();
        
        std::string sval = obj->get_rmin1();
        sval += "*"+lunit;
        double rmin1 = calc->Eval( sval );
        sval = obj->get_rmin2();
        sval += "*"+lunit;
        double rmin2 = calc->Eval( sval );
        sval = obj->get_rmax1();
        sval += "*"+lunit;
        double rmax1 = calc->Eval( sval );
        sval = obj->get_rmax2();
        sval += "*"+lunit;
        double rmax2 = calc->Eval( sval );
        sval = obj->get_z();
        sval += "*"+lunit;
        double dz = calc->Eval( sval ); dz = dz/2.;
        sval = obj->get_startphi();
        sval += "*"+aunit;
        double startphi = calc->Eval( sval );
        sval = obj->get_deltaphi();
        sval += "*"+aunit;
        double deltaphi = calc->Eval( sval );
        
//         std::cout << "rmin1:       "  << obj->get_rmin1()       << lunit
//                   << " rmin1:      "  << rmin1 << std::endl;
//         std::cout << "rmin2:       "  << obj->get_rmin2()       << lunit
//                   << " rmin2:      "  << rmin2 << std::endl;
//         std::cout << "rmax1:       "  << obj->get_rmax1()       << lunit
//                   << " rmax1:      "  << rmax1 << std::endl;
//         std::cout << "rmax2:       "  << obj->get_rmax2()       << lunit
//                   << " rmax2:      "  << rmax2 << std::endl;
//         std::cout << "z:          "  << obj->get_z()       << lunit
//                   << " dz:        "  << dz << std::endl;
//         std::cout << "startphi:   "  << obj->get_startphi()   << aunit
//                   << " startphi:  "  << startphi << std::endl;
//         std::cout << "deltaphi:   "  << obj->get_deltaphi()   << aunit
//                   << " deltaphi:  "  << deltaphi << std::endl;

        G4VSolid* newobj = new G4Cons( name, rmin1, rmax1, rmin2, rmax2, dz, startphi, deltaphi );
      
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

DECLARE_SUBSCRIBER_FACTORY(coneSubscriber)

