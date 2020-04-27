from src.system import System
from src.simulate import simulate
from src.animate import animate
from src.plot_spin_current import plot_spin_current

from src.systems.system_a import make_system_a

h_target = 0.5  # (nm)
dt = 0.5  # (fs)
electrons_per_packet = 0.00001  # (nm^-2)

system: System = make_system_a(h_target, dt, electrons_per_packet)

state_list = simulate(system, 200.0)

animate(system, state_list, 5.0)

plot_spin_current(system,  5.0, state_list)
plot_spin_current(system, 10.0, state_list)
plot_spin_current(system, 15.0, state_list)
plot_spin_current(system, 20.0, state_list)
plot_spin_current(system, 25.0, state_list)