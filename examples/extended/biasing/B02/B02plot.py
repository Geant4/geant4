# This script may be executed from a lizard session to display
# the histogram in b02.hbook. To create the histogram and the file
# b02.hbook run exampleB02 first.

tree = tf.create ("b02.hbook", "hbook",1,0)
h = tree.findH1D("10")
hplot(h)
