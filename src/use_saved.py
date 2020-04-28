import pickle

from src.system import System
from src.simulate import simulate
from src.animate import animate
from src.plot_spin_current import plot_spin_current


system = pickle.load(open("system.p", 'rb'))
state_list = pickle.load(open("state_list.p", 'rb'))

print("system made")

#animate(system, state_list, 5.0)

plot_spin_current(system,  5.0, state_list, True)
plot_spin_current(system, 10.0, state_list, True)
plot_spin_current(system, 15.0, state_list, True)
plot_spin_current(system, 20.0, state_list, True)
plot_spin_current(system, 25.0, state_list, True)
