cd dirUserRoutines
#
$FLUPRO/flutil/fff mgdraw.f
$FLUPRO/flutil/fff usrein.f
$FLUPRO/flutil/fff usreou.f
$FLUPRO/flutil/fff usrini.f
$FLUPRO/flutil/fff usrout.f
#
$FLUPRO/flutil/lfluka -m fluka mgdraw.o usrein.o usreou.o usrini.o usrout.o \
                      -o flukahp
#
mv flukahp ../.
mv flukahp.map ../.
cd ../
