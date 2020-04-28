import pickle

from src.system import System
from src.simulate import simulate
from src.animate import animate
from src.plot_spin_current import plot_spin_current

from src.systems.system_a import make_system_a
from src.systems.system_b.system_b import make_system_b

h_target = 0.1  # (nm)
dt = 1.0  # (fs)
substeps = 50
electrons_per_packet = 0.001  # (nm^-2)

system: System = make_system_b(h_target, dt, substeps, electrons_per_packet)

state_list = simulate(system, 50.0)

pickle.dump(system,     open("system.p", 'wb'))
pickle.dump(state_list, open("state_list.p", 'wb'))

animate(system, state_list, 10.0, 14.0, 16.0)
#plot_spin_current(system, 10.0, state_list, True)
