#include "commonwrapper.hh"
#include "functionwrapper.hh"

struct hazard_common* gHazard;
struct mat_common *gMat;
struct ws_common *gWs;
struct saxw_common *gSaxw;
struct light_gaus_nuc_common *gLight_gaus_nuc;
struct light_nuc_common *gLight_nuc;
struct calincl_common *gCalincl;
struct spl2_common *gSpl2;
struct dton_common *gDton;

namespace 
{
	// Helper class to initalize pointer gHazard at load time. 
	struct initializer 
	{
		static initializer _instance;
		initializer() { 
			gHazard = &hazard_; 
			gMat = &mat_;
			gWs = &ws_;
			gSaxw = &saxw_;
			gLight_gaus_nuc = &light_gaus_nuc__;
			gLight_nuc = &light_nuc__;
			gCalincl = &calincl_;
			gSpl2 = &spl2_;
			gDton = &dton_;
		}
	} _instance;
}
	
