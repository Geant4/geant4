#!/bin/sh -f
#
#  To remove test directories after a checkout.
#
find $G4INSTALL/source -name "test" -exec /bin/rm -R {} \;
#