#
#
foreach target ( testPropagateMagField testPropagateSpin testProPerpSpin testProElectroMagField )
   gmake G4TARGET=$target
end

exit

#
gmake G4TARGET=testPropagateMagField
   320  19:10   gmake G4TARGET=testPropagateSpin
   321  19:11   gmake G4TARGET=testProPerpSpin
   322  19:11   gmake G4TARGET=testProElectroMagField

