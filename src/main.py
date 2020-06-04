import pickle

from src.system import System
from src.simulate import simulate
from src.animate import animate
from src.plot_spin_current import plot_spin_current
from src.total_transfer_depthresolved import total_transfer_depthresolved

from src.systems.system_a import make_system_a
from src.systems.system_b.system_b import make_system_b
from src.systems.system_c import make_system_c

h_target = 0.25  # (nm)
dt = 0.5  # (fs)
substeps = 10
electrons_per_packet = 0.0001  # (nm^-2)

system: System = make_system_c(h_target, dt, substeps, electrons_per_packet)

state_list = simulate(system, 200.0)

pickle.dump(system,     open("system.p", 'wb'))
pickle.dump(state_list, open("state_list.p", 'wb'))

animate(system, state_list, 10.0, 0.0, 40.0)

#plot_spin_current(system, state_list, 16.0, False)
#total_transfer_depthresolved(system, state_list, 15.0, 30.0, 0.01)