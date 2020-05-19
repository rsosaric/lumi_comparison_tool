import pickle as pckl
import matplotlib.pyplot as plt
import optparse

p = optparse.OptionParser()
(options, args) = p.parse_args()

input_path = args[0]

fig = pckl.load(open(input_path,  'rb'))
plt.show()
