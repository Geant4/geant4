#include <stdio.h>
#include <fenv.h>
void fpu_ ()
{
  fesetenv (FE_NOMASK_ENV);
}
