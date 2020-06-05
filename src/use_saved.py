import matplotlib.pyplot as plt
import pickle

from src.system import System
from src.simulate import simulate
from src.animate import animate
from src.plot_spin_current import plot_spin_current
from src.injected import injected
from src.plot_fourier import plot_fourier
from src.total_transfer_depthresolved import total_transfer_depthresolved
from src.plot_gamma import plot_gamma

system = pickle.load(open("system.p", 'rb'))
state_list = pickle.load(open("state_list.p", 'rb'))

print("system made")

#animate(system, state_list, 10.0, 0.0, 200.0)

plot_gamma(system, state_list, 150.0, 0.0, 40.0, -0.8, 0.2)

plot_spin_current(system, state_list, 40.0, False)
plot_spin_current(system, state_list, 50.0, False)
plot_spin_current(system, state_list, 60.0, False)
plot_spin_current(system, state_list, 70.0, False)
plot_spin_current(system, state_list, 80.0, False)
plot_spin_current(system, state_list, 90.0, False)

#r, theta = injected(system, state_list, 20.0, 0.01)
#print(r)
#print(theta)

#z = 17.0
#plot_spin_current(system, state_list, z, True)
#plot_fourier(system, state_list, z, 0.1)

#total_transfer_depthresolved(system, state_list, 15.0, 40.0, 0.01)
