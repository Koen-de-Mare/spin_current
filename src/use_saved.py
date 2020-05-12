import matplotlib.pyplot as plt
import pickle

from src.system import System
from src.simulate import simulate
from src.animate import animate
from src.plot_spin_current import plot_spin_current
from src.injected import injected
from src.plot_fourier import plot_fourier
from src.total_transfer_depthresolved import total_transfer_depthresolved

system = pickle.load(open("system.p", 'rb'))
state_list = pickle.load(open("state_list.p", 'rb'))

print("system made")

#animate(system, state_list, 10.0, 0.0, 100.0)

#r, theta = injected(system, state_list, 20.0, 0.01)
#print(r)
#print(theta)

#z = 30.0
#plot_spin_current(system, state_list, z, False)
#plot_fourier(system, state_list, z, 0.1)

total_transfer_depthresolved(system, state_list, 15.0, 30.0, 0.01)
