#include "CLHEP/Random/DoubConv.h"

#include "CLHEP/Random/RandExpZiggurat.h"

#include <iostream>
#include <cmath>	// for std::log()

namespace CLHEP {

CLHEP_THREAD_LOCAL unsigned long RandExpZiggurat::kn[128], RandExpZiggurat::ke[256];
CLHEP_THREAD_LOCAL float RandExpZiggurat::wn[128],RandExpZiggurat::fn[128],RandExpZiggurat::we[256],RandExpZiggurat::fe[256];
CLHEP_THREAD_LOCAL bool RandExpZiggurat::ziggurat_is_init = false;

std::string RandExpZiggurat::name() const {return "RandExpZiggurat";}

HepRandomEngine & RandExpZiggurat::engine() {return *localEngine;}

RandExpZiggurat::~RandExpZiggurat() {
}

RandExpZiggurat::RandExpZiggurat(const RandExpZiggurat& right) : HepRandom(right),defaultMean(right.defaultMean)
{
}

double RandExpZiggurat::operator()()
{
  return fire( defaultMean );
}

void RandExpZiggurat::shootArray( const int size, float* vect, float mean )
{
   for (int i=0; i<size; ++i) vect[i] = shoot(mean);
}

void RandExpZiggurat::shootArray( const int size, double* vect, double mean )
{
   for (int i=0; i<size; ++i) vect[i] = shoot(mean);
}

void RandExpZiggurat::shootArray(HepRandomEngine* anEngine, const int size, float* vect, float mean )
{
   for (int i=0; i<size; ++i) vect[i] = shoot(anEngine, mean);
}

void RandExpZiggurat::shootArray(HepRandomEngine* anEngine, const int size, double* vect, double mean )
{
   for (int i=0; i<size; ++i) vect[i] = shoot(anEngine, mean);
}

void RandExpZiggurat::fireArray( const int size, float* vect)
{
   for (int i=0; i<size; ++i) vect[i] = fire( defaultMean );
}

void RandExpZiggurat::fireArray( const int size, double* vect)
{
   for (int i=0; i<size; ++i) vect[i] = fire( defaultMean );
}

void RandExpZiggurat::fireArray( const int size, float* vect, float mean )
{
   for (int i=0; i<size; ++i) vect[i] = fire( mean );
}

void RandExpZiggurat::fireArray( const int size, double* vect, double mean )
{
   for (int i=0; i<size; ++i) vect[i] = fire( mean );
}

std::ostream & RandExpZiggurat::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  std::vector<unsigned long> t(2);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  t = DoubConv::dto2longs(defaultMean);
  os << defaultMean << " " << t[0] << " " << t[1] << "\n";
  os.precision(pr);
  return os;
#ifdef REMOVED
  int pr=os.precision(20);
  os << " " << name() << "\n";
  os << defaultMean << "\n";
  os.precision(pr);
  return os;
#endif
}

std::istream & RandExpZiggurat::get ( std::istream & is ) {
  std::string inName;
  is >> inName;
  if (inName != name()) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Mismatch when expecting to read state of a "
    	      << name() << " distribution\n"
	      << "Name found was " << inName
	      << "\nistream is left in the badbit state\n";
    return is;
  }
  if (possibleKeywordInput(is, "Uvec", defaultMean)) {
    std::vector<unsigned long> t(2);
    is >> defaultMean >> t[0] >> t[1]; defaultMean = DoubConv::longs2double(t); 
    return is;
  }
  // is >> defaultMean encompassed by possibleKeywordInput
  return is;
}


float RandExpZiggurat::ziggurat_efix(unsigned long jz,HepRandomEngine* anEngine)
{ 
  if(!ziggurat_is_init) ziggurat_init();

  unsigned long iz=jz&255;
  
  float x;
  for(;;)
  {  
    if(iz==0) return (7.69711-std::log(ziggurat_UNI(anEngine)));          /* iz==0 */
    x=jz*we[iz]; 
    if( fe[iz]+ziggurat_UNI(anEngine)*(fe[iz-1]-fe[iz]) < std::exp(-x) ) return (x);

    /* initiate, try to exit for(;;) loop */
    jz=ziggurat_SHR3(anEngine);
    iz=(jz&255);
    if(jz<ke[iz]) return (jz*we[iz]);
  }
}

bool RandExpZiggurat::ziggurat_init()
{  
  const double rzm1 = 2147483648.0, rzm2 = 4294967296.;
  double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
  double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
  int i;

/* Set up tables for RNOR */
  q=vn/std::exp(-.5*dn*dn);
  kn[0]=(unsigned long)((dn/q)*rzm1);
  kn[1]=0;

  wn[0]=q/rzm1;
  wn[127]=dn/rzm1;

  fn[0]=1.;
  fn[127]=std::exp(-.5*dn*dn);

  for(i=126;i>=1;i--) {
    dn=std::sqrt(-2.*std::log(vn/dn+std::exp(-.5*dn*dn)));
    kn[i+1]=(unsigned long)((dn/tn)*rzm1);
    tn=dn;
    fn[i]=std::exp(-.5*dn*dn);
    wn[i]=dn/rzm1;
  }

/* Set up tables for REXP */
  q = ve/std::exp(-de);
  ke[0]=(unsigned long)((de/q)*rzm2);
  ke[1]=0;

  we[0]=q/rzm2;
  we[255]=de/rzm2;

  fe[0]=1.;
  fe[255]=std::exp(-de);

  for(i=254;i>=1;i--) {
    de=-std::log(ve/de+std::exp(-de));
    ke[i+1]= (unsigned long)((de/te)*rzm2);
    te=de;
    fe[i]=std::exp(-de);
    we[i]=de/rzm2;
  }
  ziggurat_is_init=true;
  return true;
}

}  // namespace CLHEP
