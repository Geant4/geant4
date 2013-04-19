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
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02SteppingAction.hh"
#include "G4SteppingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleTypeFilter : public G4VLoggingFilter {
private:
    G4ParticleDefinition* part;
public:
    G4ParticleTypeFilter( G4ParticleDefinition* pd ) : part(pd) { }
    G4String Name() const { return "ParticleTypeFilter_"+ (part->GetParticleName()); }
    G4ParticleTypeFilter* Clone() const { return new G4ParticleTypeFilter( part ); }
    G4bool AcceptStep( const G4Step* aStep ) {
        if ( aStep->GetTrack()->GetParticleDefinition() == part )
            return true;
        else
            return false;
    }
};

#include <sstream>
#include "G4SystemOfUnits.hh"
class G4RangeFilter : public G4VLoggingFilter {
protected:
  G4double min;
  G4double max;
public:
    G4RangeFilter( G4double aMin , G4double aMax ) : min(aMin), max(aMax) { }
    G4bool CheckQuantity( G4double q )
    {
        if ( q >= min && q < max ) return true;
        else return false;
    }
};

class G4EnergyRangeFilter : public G4RangeFilter {
public:
    G4EnergyRangeFilter( G4double minInMeV , G4double maxInMeV ) : G4RangeFilter(minInMeV,maxInMeV) {}
    G4String Name() const {
        std::ostringstream msg;
        msg << "EnergyRangeFilter_"<< min <<"-"<< max <<"_MeV";
        return msg.str();
    }
    G4EnergyRangeFilter* Clone() const { return new G4EnergyRangeFilter(min,max); }
    G4bool AcceptStep( const G4Step* aStep )
    {
        G4double ene = aStep->GetPostStepPoint()->GetKineticEnergy();
        return CheckQuantity(ene/MeV);
    }
};

class G4StepLengthRangeFilter : public G4RangeFilter {
public:
    G4StepLengthRangeFilter( G4double min_mm , G4double max_mm ) : G4RangeFilter(min_mm,max_mm) {}
    G4String Name() const {
        std::ostringstream msg;
        msg << "StepLengthRangeFilter_"<< min <<"-"<< max <<"_mm";
        return msg.str();
    }
    G4StepLengthRangeFilter* Clone() const { return new G4StepLengthRangeFilter(min,max); }
    G4bool AcceptStep( const G4Step* aStep )
    {
        G4double l = aStep->GetStepLength();
        return CheckQuantity(l/mm);
    }
};

class G4MaterialFilter : public G4VLoggingFilter {
private:
    const G4String mat;
public:
    G4MaterialFilter( const G4String& matName ) : mat(matName) { }
    G4String Name() const { return "MaterialFilter-"+mat; }
    G4MaterialFilter* Clone() const { return new G4MaterialFilter(mat); }
    G4bool AcceptStep( const G4Step* aStep )
    {
        G4StepPoint* pre = aStep->GetPreStepPoint();
        G4Material*  stepmat = pre->GetMaterial();
        if ( stepmat->GetName() == mat ) return true;
        else return false;
    }
};

class G4VolumeFilter : public G4VLoggingFilter {
private:
    const G4String volname;
public:
    G4VolumeFilter( const G4String& volName ) : volname(volName ) { }
    G4String Name() const { return "VolumeFilter-"+volname; }
    G4VolumeFilter* Clone() const { return new G4VolumeFilter(volname); }
    G4bool AcceptStep(const G4Step* aStep )
    {
        const G4String& n = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
        if ( n == volname ) return true;
        else return false;
    }
};

class G4ProcessFilter : public G4VLoggingFilter {
public:
    virtual G4bool AcceptProc( const G4VProcess* aProc) = 0 ;
    G4bool AcceptStep(const G4Step* aStep )
    {
        const G4VProcess* proc = aStep->GetPostStepPoint()->GetProcessDefinedStep();
        return AcceptProc(proc);
    }
};

class G4ProcessTypeFilter : public G4ProcessFilter {
private:
    G4int proctype;
    G4int procsubtype;
public:
    G4ProcessTypeFilter( G4int procType , G4int procSubType = -1 ) : proctype(procType) , procsubtype(procSubType) { }
    G4String Name() const
    {
        std::ostringstream msg;
        msg<<"ProcessTypeFilter_"<<proctype;
        if ( procsubtype > -1 ) msg<<"_"<<procsubtype;
        return msg.str();
    }
    G4bool AcceptProc( const G4VProcess* proc )
    {
        if ( proctype == proc->GetProcessType() )
        {
            if ( procsubtype > -1 ) {
                if ( procsubtype == proc->GetProcessSubType() ) { return true; }
                else { return false; }
            } else { return true; }
        } else { return false; }
    }
    G4ProcessTypeFilter* Clone() const { return new G4ProcessTypeFilter(proctype,procsubtype); }
};

class G4ProcessNameFilter : public G4ProcessFilter {
private:
    G4String pname;
public:
    G4ProcessNameFilter(const G4String& pn ) : pname(pn) { }
    G4String Name() const { return "ProcessNameFilter_"+pname; }
    G4bool AcceptProc(const G4VProcess* proc )
    {
        if ( pname == proc->GetProcessName() ) return true;
        else return false;
    }
};

G4DistributionInfo::~G4DistributionInfo()
{
    while (!infos.empty())
    {
        G4StepInfo* e = *(infos.begin());
        infos.pop_front();
        delete e;
    }
}



G4bool G4DistributionInfo::LogStep(const G4Step* aStep, G4double quantity)
{
    ++counter;
    G4bool result = true;
    for ( std::list<G4StepInfo*>::const_iterator it = infos.begin() ;
         it != infos.end() ;
        ++it )
    {
        result &= (*it)->LogStep(aStep,quantity);
    }
    return result;
}

G4DistributionInfo * dist = 0;
#include "G4MuonPlus.hh"
#include "G4Positron.hh"
ExN02SteppingAction::ExN02SteppingAction()
{
    dist = new G4DistributionInfo("Ene");
    //Energy distribution of e+
    G4StepInfo * info = new G4StepInfo("EneDist1");
    info->AddFilter(new G4ParticleTypeFilter(G4Positron::Definition() ) );
    dist->AddInfoLogger(info);
    //Energy distribution of parciles with 100MeV<=Ekin<200MeV
    info = new G4StepInfo("EneDist2");
    info->AddFilter( new G4EnergyRangeFilter(100,200) );
    dist->AddInfoLogger(info);
    //Energy distribution of steps with 0.9mm<=L<1mm
    info = new G4StepInfo("EneDist3");
    info->AddFilter(new G4StepLengthRangeFilter(0.9,1) );
    dist->AddInfoLogger(info);
    //Combine: like previous and for energy <100
    G4StepInfo* info2 = new G4StepInfo(*info);
    info2->SetName("EneDistComb");
    info2->AddFilter(new G4EnergyRangeFilter(0,100) );
    dist->AddInfoLogger(info2);
    //Material
    info = new G4StepInfo("EneDist4");
    info->AddFilter(new G4MaterialFilter("Lead") );
    dist->AddInfoLogger(info);
    //Volume
    info = new G4StepInfo("EneDist5");
    info->AddFilter(new G4VolumeFilter("Chamber"));
    dist->AddInfoLogger(info);
    //Process type
    info = new G4StepInfo("EneDist6");
    info->AddFilter(new G4ProcessTypeFilter(12) );
    dist->AddInfoLogger(info);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    /*G4bool result =*/ dist->LogStep(aStep,aStep->GetTotalEnergyDeposit());
}

#include <fstream>
void ExN02SteppingAction::Dump() const
{
    G4cout<<"Writing monitoring stuff in data.xml"<<G4endl;
    std::ofstream fout;
    fout.open("data.xml");
    fout<<"<statistics>\n";
    fout<<*dist;
    fout<<"</statistics>";
    fout.close();
}

ExN02SteppingAction::~ExN02SteppingAction()
{
    delete dist;
    dist=0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//#define STEPINFODEBUG
#ifdef STEPINFODEBUG
G4bool G4StepInfo::debug = true;
#else
G4bool G4StepInfo::debug = false;
#endif


std::ostream& operator<<(std::ostream& out, const G4StepInfo& stepInfo )
{
    stepInfo.Dump(out);
    return out;
}

std::ostream& operator<<(std::ostream& out, const G4VLoggingFilter& logfileter)
{
    logfileter.Dump(out);
    return out;
}

std::ostream& operator<<(std::ostream& out,const G4DistributionInfo& mydist )
{
    mydist.Dump(out);
    return out;
}

void G4VLoggingFilter::Dump(std::ostream& out) const
{
    out<<"<G4VLoggingFilter name=\""<<Name()<<"\"/>\n";
}

void G4StepInfo::Clear()
{
    vals.clear();
    while (! filters.empty() )
    {
        G4VLoggingFilter* f = *(filters.begin());
        filters.pop_front();
        delete f;
    }
}

void G4StepInfo::Dump(std::ostream& out) const
{
    out<<"<G4StepInfo name=\""<<name<<"\">\n";
    out<<"<results>\n<counter>"<<counter<<"</counter>\n";
    out<<"<mean>"<<Mean()<<"</mean>\n"<<"<variance>"<<Variance()<<"</variance>\n"<<"<excesskurtosis>"<<ExcKurtosis()<<"</excesskurtosis>\n"<<"<skewness>"<<Skewness()<<"</skewness>\n";
    out<<"</results>\n";
    out<<"<filters>\n";
    for ( std::list<G4VLoggingFilter*>::const_iterator it = filters.begin() ;
         it != filters.end() ; ++it )
        out<<*(*it);
    out<<"</filters>\n";
    if ( debug )
    {
        out<<"<data>";
        for ( std::list<float>::const_iterator it = vals.begin();
             it != vals.end() ; ++it )
            out<<*it<<" ";
        out<<"\n</data/>";
    }
    out<<"</G4StepInfo>\n";
}

void G4DistributionInfo::Dump(std::ostream& out) const
{
    out<<"<G4DistributionInfo name=\""<<name<<"\">\n";
    out<<"<counter>"<<counter<<"</counter>\n";
    out<<"<loggers>\n";
    for ( std::list<G4StepInfo*>::const_iterator it = infos.begin() ;
         it != infos.end() ;
         ++it ) out<<*(*it);
    out<<"</loggers>\n";
    out<<"</G4DistributionInfo>\n";
}


G4StepInfo::G4StepInfo(const char* n)
    : name(n),counter(0),sum(0),sum2(0),sum3(0),sum4(0)
{}

G4StepInfo::G4StepInfo(const G4StepInfo& rhs )
{
    if ( this == &rhs ) return;
    Clear();
    name = rhs.name;
    vals = rhs.vals;
    sum = rhs.sum;
    sum2 = rhs.sum2;
    sum3 = rhs.sum3;
    sum4 = rhs.sum4;
    counter = rhs.counter;
    for ( std::list<G4VLoggingFilter*>::const_iterator it = rhs.filters.begin() ;
         it != rhs.filters.end() ; ++it )
    {
        filters.push_back( (*it)->Clone() );
    }
    for ( std::list<float>::const_iterator it = rhs.vals.begin() ;
         it != rhs.vals.end() ; ++it )
    {
        vals.push_back(*it);
    }
}

void G4StepInfo::AddFilter(G4VLoggingFilter* f) { filters.push_back(f); }

G4StepInfo::~G4StepInfo()
{
    if (debug )
    {
        std::ostringstream msg;
        Dump(msg);
        G4cout<<msg<<G4endl;
    }
    Clear();
}

G4bool G4StepInfo::LogStep( const G4Step* aStep , G4double value )
{
    G4bool result = true;
    for ( std::list<G4VLoggingFilter*>::iterator it = filters.begin() ;
         it!=filters.end() ; ++it)
    {
        result &= (*it)->AcceptStep(aStep);
        //no break if false, maybe filters do some logging on their own
    }
    if ( result )
    {
        ++counter;
        sum += value;
        sum2 += (value*value);
        sum3 += (value*value*value);
        sum4 += (value*value*value*value);
        if ( debug ) vals.push_back(value);
    }
    return result;
}

G4double G4StepInfo::Mean() const { return G4double(sum)/G4double(counter); }

G4double G4StepInfo::Variance() const
{
    if (counter<2) return -DBL_MAX;
    G4double mean = Mean();
    G4double mean2 = mean*mean;
    G4double s1 = sum2 + counter*mean2 - 2*mean*sum;
    return s1/(counter-1);
}

G4double G4StepInfo::ExcKurtosis() const
{
    if (counter<4) return -DBL_MAX;
    G4double counter2 = counter*counter;
    G4double mom4= -3*sum*sum*sum*sum/(counter2*counter2) \
                   +6*sum*sum*sum2/(counter2*counter) \
                    -4*sum*sum3/counter2 + sum4/counter;
    G4double var = sum2/counter-(sum*sum)/(counter2);//Sample variance, not population
    if ( var < 1E-10 ) return -DBL_MAX;
    G4double corr = 3*(counter-1)*(counter-1);
    corr /= ( (counter-2)*(counter-3));
    G4double fact = (counter-1)*(counter+1)/( (counter-2)*(counter-3) );
    return fact*mom4/(var*var)-corr;
}

#include <cmath>
G4double G4StepInfo::Skewness() const
{
    if (counter<3) return -DBL_MAX;
    G4double counter2 = counter*counter;
    G4double coeff = sqrt( counter*(counter-1) )/(counter-2);
    G4double mom3 = 2*sum*sum*sum/(counter*counter2) \
                    -3*sum*sum2/counter2 + sum3/counter;
    G4double var = sum2/counter-(sum*sum)/(counter2);//Sample variance, not population
    if ( var < 1E-10 ) return -DBL_MAX;
    return coeff*mom3/pow(var,1.5);
}
