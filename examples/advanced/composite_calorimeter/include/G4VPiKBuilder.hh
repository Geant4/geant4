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
#ifndef G4VPiKBuilder_h
#define G4VPiKBuilder_h

#include "G4HadronElasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"

class G4VPiKBuilder
{
  public:
    G4VPiKBuilder() {}
    virtual ~G4VPiKBuilder() {}
    virtual void Build(G4HadronElasticProcess & aP) = 0;
    virtual void Build(G4PionPlusInelasticProcess & aP) = 0;
    virtual void Build(G4PionMinusInelasticProcess & aP) = 0;
    virtual void Build(G4KaonPlusInelasticProcess & aP) = 0;
    virtual void Build(G4KaonMinusInelasticProcess & aP) = 0;
    virtual void Build(G4KaonZeroLInelasticProcess & aP) = 0;
    virtual void Build(G4KaonZeroSInelasticProcess & aP) = 0;
};
// 2002 by J.P. Wellisch

#endif
