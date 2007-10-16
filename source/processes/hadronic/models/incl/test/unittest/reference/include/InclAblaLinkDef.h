
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// Common blocks
#pragma link C++ struct hazard_common;
#pragma link C++ struct mat_common;
#pragma link C++ struct ws_common;
#pragma link C++ struct saxw_common;
#pragma link C++ struct light_gaus_nuc_common;
#pragma link C++ struct light_nuc_common;
#pragma link C++ struct calincl_common;
#pragma link C++ struct spl2_common;
#pragma link C++ struct dton_common;

#pragma link C++ global gHazard;
#pragma link C++ global gMat;
#pragma link C++ global gWs;
#pragma link C++ global gSaxw;
#pragma link C++ global gLight_gaus_nuc;
#pragma link C++ global gLight_nuc;
#pragma link C++ global gCalincl;
#pragma link C++ global gSpl2;
#pragma link C++ global gDton;

// Subroutines and functions
// init_incl42:
#pragma link C++ function init_incl__;
#pragma link C++ function init_mat__;
#pragma link C++ function integ_;
#pragma link C++ function flin2_;
#pragma link C++ function dens_deut__;
#pragma link C++ function spl2ab_;
#pragma link C++ function dens_;

#pragma link C++ function wsax_;
#pragma link C++ function derivwsax_;
#pragma link C++ function dmho_;
#pragma link C++ function derivmho_;
#pragma link C++ function derivgaus_;
#pragma link C++ function flin_;
#pragma link C++ function flin2_;
#pragma link C++ function deutv_;
#pragma link C++ function fm2_;
#pragma link C++ function dens_;
#pragma link C++ function splineab_;

//incl4.2:
#pragma link C++ function ribm_;
#pragma link C++ function rgauss_;
#pragma link C++ function texp_;
#pragma link C++ function sig_reac__;
#pragma link C++ function radi_us__;
#pragma link C++ function xabs2_;
#pragma link C++ function coulomb_transm__;
#pragma link C++ function clmb1_;
#pragma link C++ function clmb2_;
#pragma link C++ function force_abs__;

//fission_distri:

//abla:
#pragma link C++ function init_evapora__;
#pragma link C++ function lpoly_;
#pragma link C++ function barfit_;
#pragma link C++ function bipol_;
#pragma link C++ function cram_;
#pragma link C++ function spdef_;
// CERNLIB:
#pragma link C++ function ranf_;

#endif
