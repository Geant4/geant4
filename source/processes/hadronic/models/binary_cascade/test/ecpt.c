

#include <fenv.h>

 
 int ecpt()
{
     
//     fesetenv (FE_NOMASK_ENV);
//     feclearexcept(FE_INEXACT);

     feenableexcept(FE_INVALID);
     feenableexcept(FE_DIVBYZERO);
     return 0;
}
