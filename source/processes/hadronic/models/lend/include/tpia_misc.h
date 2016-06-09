/*
# <<BEGIN-copyright>>
# Copyright (c) 2010, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory 
# Written by Bret R. Beck, beck6@llnl.gov. 
# CODE-461393
# All rights reserved. 
#  
# This file is part of GIDI. For details, see nuclear.llnl.gov. 
# Please also read the "Additional BSD Notice" at nuclear.llnl.gov. 
# 
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met: 
#
#      1) Redistributions of source code must retain the above copyright notice, 
#         this list of conditions and the disclaimer below.
#      2) Redistributions in binary form must reproduce the above copyright notice, 
#         this list of conditions and the disclaimer (as noted below) in the 
#          documentation and/or other materials provided with the distribution.
#      3) Neither the name of the LLNS/LLNL nor the names of its contributors may be 
#         used to endorse or promote products derived from this software without 
#         specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT 
# SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS 
# OR SERVICES;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
# AND ON  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
# <<END-copyright>>
*/

#if defined __cplusplus
    namespace GIDI {
    extern "C" {
#endif

#define tpia_misc_pointerToAttributeIfAllOk2( smr, element, required, attributes, name ) \
    tpia_misc_pointerToAttributeIfAllOk( smr, element, NULL, required, attributes, name, __FILE__, __LINE__ )
#define tpia_misc_pointerToAttributeIfAllOk3( smr, path, required, attributes, name ) \
    tpia_misc_pointerToAttributeIfAllOk( smr, NULL, path, required, attributes, name, __FILE__, __LINE__ )

char *tpia_misc_pointerToAttributeIfAllOk( statusMessageReporting *smr, xData_element *element, const char *path, int required, 
        xData_attributionList *attributes, const char *name, const char *file, int line );
int tpia_misc_setMessageError_Element( statusMessageReporting *smr, void *userInterface, xData_element *element, const char *file, int line, int code, 
    const char *fmt, ... );
xData_Int tpia_misc_binarySearch( xData_Int n, double *ds, double d );

tpia_EqualProbableBinSpectrum *tpia_misc_getEqualProbableBin( statusMessageReporting *smr, xData_element *parent, xData_Int *n, xData_Int *nBins );
tpia_EqualProbableBinSpectrum *tpia_misc_getEqualProbableBins( statusMessageReporting *smr, xData_element *parent, const char *name, xData_Int nBins, 
    xData_Int *n );

#if defined __cplusplus
    }
    }
#endif
