#!/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/PublicDomainPackages/2.0.0/bin/python2.2
# #!/usr/bin/env python2


#  This example needs python2.2 http://www.python.org or later.
#  
#  This script may be used without or with ANAPHE e.g. lizard.
#  Without lizrd:  just execute this file.
#  With lizrd:     start lizard and than  do
#                    execfile("B03RunApplication.py")

import string


execfile("B03Application.py")
#executing the application could be done "by hand"

Bapp = B03Application(15)
Bapp.initializeApplication()
Bapp.setupParallelGeometry()
Bapp.createIstoreAndScorer()
Bapp.setupSampler()

Bapp.run(15)

# using lizard and ANAPHE to create a histogram
# check if lizrd is loaded by looking for the af instance of the AnalysisFactory
if 'af' in  dir():     
    tr = tf.create()
    hf = af.createHistogramFactory(tr)
    n = Bapp.nevents
    minm = int(min(Bapp.messurements))-0.5
    maxm = int(max(Bapp.messurements)+1)+0.5
    h = hf.createHistogram1D ("10","Tracks entering in " +  
                              '%(n)d' % vars() +  
                              " events",30, minm, maxm)
    for v in Bapp.messurements:
        h.fill(v)

    hplot(h)

