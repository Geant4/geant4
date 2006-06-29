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
//
// $Id: SolidAnalyser.cc,v 1.6 2006-06-29 17:23:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// SolidAnalyser
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#include "SolidAnalyser.hh"

#include <sstream>

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"

SolidAnalyser * SolidAnalyser::theInstance = 0;

SolidAnalyser::SolidAnalyser()
{
}


SolidAnalyser::~SolidAnalyser()
{
   delete theInstance;
}


SolidAnalyser * SolidAnalyser::GetSolidAnalyser()
{
   if (! theInstance) 
     theInstance = new SolidAnalyser();
   
   return theInstance;  
}


G4int SolidAnalyser::GetParam(const G4VSolid * s,
             std::vector<std::pair<G4String,G4double> > & result ) const
{
   const G4Box * box = dynamic_cast<const G4Box*>(s);   
   if ( box ) return GetParam(box,result);
   
   const G4Cons * cons = dynamic_cast<const G4Cons*>(s);  
   if ( cons ) return GetParam(cons,result);
   
   const G4Polycone * pcon = dynamic_cast<const G4Polycone*>(s);  
   if ( pcon ) return GetParam(pcon,result);
   
   const G4Polyhedra * phed = dynamic_cast<const G4Polyhedra*>(s);  
   if ( phed ) return GetParam(phed,result);
   
   const G4Trap * trap = dynamic_cast<const G4Trap*>(s);  
   if ( trap ) return GetParam(trap,result);
     
   const G4Trd * trd = dynamic_cast<const G4Trd*>(s);  
   if ( trd ) return GetParam(trd,result);
     
   const G4Tubs * tubs = dynamic_cast<const G4Tubs*>(s);
   if ( tubs ) return GetParam(tubs,result);
   
   return NotImplemented(s,result);  
}

// default parameters for not implemented solids
G4int SolidAnalyser::NotImplemented(const G4VSolid * s,
         std::vector<std::pair<G4String,G4double> > & result) const
{
   G4String temp = s->GetEntityType() + G4String(": not impl.");
   result.push_back(std::make_pair(temp,-1.0));
   return -1;
}

// G4Box
G4int SolidAnalyser::GetParam(const G4Box * s,
          std::vector<std::pair<G4String,G4double> > & result) const
{
   result.push_back(std::make_pair(G4String("x/2"), s->GetXHalfLength()));
   result.push_back(std::make_pair(G4String("y/2"), s->GetYHalfLength()));
   result.push_back(std::make_pair(G4String("z/2"), s->GetZHalfLength()));
   return 3;
}

// G4Cons
G4int SolidAnalyser::GetParam(const G4Cons * s,
           std::vector<std::pair<G4String,G4double> > & result) const
{
   result.push_back(std::make_pair(G4String("rInZ-"),
                    s->GetInnerRadiusMinusZ()));
   result.push_back(std::make_pair(G4String("rInZ+"),
                    s->GetInnerRadiusPlusZ()));
   result.push_back(std::make_pair(G4String("rOutZ-"),
                    s->GetOuterRadiusMinusZ()));
   result.push_back(std::make_pair(G4String("rOutZ+"),
                    s->GetOuterRadiusPlusZ()));
   result.push_back(std::make_pair(G4String("z/2"),
                    s->GetZHalfLength()));
   result.push_back(std::make_pair(G4String("startPhi"),
                    s->GetStartPhiAngle()));
   result.push_back(std::make_pair(G4String("deltaPhi"),
                    s->GetDeltaPhiAngle()));
   return 7;
}

// G4Polycone
G4int SolidAnalyser::GetParam(const G4Polycone * s,
           std::vector<std::pair<G4String,G4double> > & result) const
{
   result.push_back(std::make_pair(G4String("startPhi"), s->GetStartPhi()));
   result.push_back(std::make_pair(G4String("endPhi"), s->GetEndPhi()));
   
   G4int nr = s->GetNumRZCorner();   
   G4String temp("nrRZ");
   result.push_back(std::make_pair(temp, G4double(nr)));
   for (G4int i=0; i<nr; i++)
   {
      std::ostringstream sstr_r, sstr_z;      
      G4PolyconeSideRZ sideRz = s->GetCorner(i);
      sstr_z << "z_" << i << '\0';
      sstr_r << "r_" << i << '\0';
      G4String z_str(sstr_z.str());
      G4String r_str(sstr_r.str());
      result.push_back(std::make_pair(z_str, sideRz.z));
      result.push_back(std::make_pair(r_str, sideRz.r));
   }
   return 3 + 2*nr;
}

// G4Polyhedra
G4int SolidAnalyser::GetParam(const G4Polyhedra * s,
         std::vector<std::pair<G4String,G4double> > & result) const
{
   result.push_back(std::make_pair(G4String("startPhi"), s->GetStartPhi()));
   result.push_back(std::make_pair(G4String("endPhi"), s->GetEndPhi()));
   result.push_back(std::make_pair(G4String("nrSide"),
                    G4double(s->GetNumSide())));
   
   G4int nr = s->GetNumRZCorner();   

   result.push_back(std::make_pair(G4String("nrRZ"), G4double(nr)));
      
   for (G4int i=0; i<nr; i++)
   {
      std::ostringstream sstr_r, sstr_z;      
      G4PolyhedraSideRZ sideRz = s->GetCorner(i);
      sstr_z << "z_" << i << '\0';
      sstr_r << "r_" << i << '\0';
      G4String z_str(sstr_z.str());
      G4String r_str(sstr_r.str());
      result.push_back(std::make_pair(z_str, sideRz.z));
      result.push_back(std::make_pair(r_str, sideRz.r));
   }   
   return 4 + 2*nr;
}

// G4Trap
G4int SolidAnalyser::GetParam(const G4Trap * s,
            std::vector<std::pair<G4String,G4double> > & result) const
{
   result.push_back(std::make_pair(G4String("z/2"), s->GetZHalfLength()));

   result.push_back(std::make_pair(G4String("y/2_1"), s->GetYHalfLength1()));
   result.push_back(std::make_pair(G4String("x/2_1"), s->GetXHalfLength1()));
   result.push_back(std::make_pair(G4String("x/2_2"), s->GetXHalfLength2()));
   result.push_back(std::make_pair(G4String("tanAlpha_1"),s->GetTanAlpha1()));
   
   result.push_back(std::make_pair(G4String("y/2_2"), s->GetYHalfLength2()));
   result.push_back(std::make_pair(G4String("x/2_3"), s->GetXHalfLength3()));
   result.push_back(std::make_pair(G4String("x/2_4"), s->GetXHalfLength4()));
   result.push_back(std::make_pair(G4String("tanAlpha_2"),s->GetTanAlpha2()));

   return 9;
}

// G4Trd
G4int SolidAnalyser::GetParam(const G4Trd * s,
             std::vector<std::pair<G4String,G4double> > & result) const
{
   result.push_back(std::make_pair(G4String("x/2_1"), s->GetXHalfLength1()));
   result.push_back(std::make_pair(G4String("x/2_2"), s->GetXHalfLength2()));
   result.push_back(std::make_pair(G4String("y/2_1"), s->GetYHalfLength1()));
   result.push_back(std::make_pair(G4String("y/2_2"), s->GetYHalfLength2()));
   result.push_back(std::make_pair(G4String("z/2"), s->GetZHalfLength()));

   return 5;
}

// G4Tubs
G4int SolidAnalyser::GetParam(const G4Tubs * s,
           std::vector<std::pair<G4String,G4double> > & result) const
{
   result.push_back(std::make_pair(G4String("rIn"), s->GetInnerRadius()));
   result.push_back(std::make_pair(G4String("rOut"), s->GetOuterRadius()));
   result.push_back(std::make_pair(G4String("z/2"), s->GetZHalfLength()));
   result.push_back(std::make_pair(G4String("startPhi"),
                    s->GetStartPhiAngle()));
   result.push_back(std::make_pair(G4String("deltaPhi"),
                    s->GetDeltaPhiAngle()));
   
   return 5;
}

std::ostream & operator<<(std::ostream& flux, 
                            std::vector<std::pair<G4String,G4double> > & v)
{
   std::vector<std::pair<G4String,G4double> >::iterator it = v.begin();
   while ( it != v.end() )
   {
      flux << (*it).first << '=' << (*it).second << "<br/>" << G4endl;
      it++;
   }   
   return flux;
}
