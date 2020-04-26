from src.system import System
from src.animate import animate
from src.simulate import simulate

from src.systems.system_a import make_system_a

h_target = 0.5  # (nm)
dt = 0.5  # (fs)
electrons_per_packet = 0.001  # (nm^-2)

system: System = make_system_a(h, dt, electrons_per_packet)

state_list = simulate(system, 50.0)

animate(system, state_list, 5.0)