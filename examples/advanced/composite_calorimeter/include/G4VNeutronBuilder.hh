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
#ifndef G4VNeutronBuilder_h
#define G4VNeutronBuilder_h

class G4HadronElasticProcess;
class G4HadronFissionProcess;
class G4HadronCaptureProcess;
class G4NeutronInelasticProcess;

class G4VNeutronBuilder
{
  public:
    G4VNeutronBuilder() {}
    virtual ~G4VNeutronBuilder() {}
    virtual void Build(G4HadronElasticProcess & aP) = 0;
    virtual void Build(G4HadronFissionProcess & aP) = 0;
    virtual void Build(G4HadronCaptureProcess & aP) = 0;
    virtual void Build(G4NeutronInelasticProcess & aP) = 0;
};
// 2002 by J.P. Wellisch

#endif
