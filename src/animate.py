import matplotlib.pyplot as plt
import matplotlib.animation as animation

from src.system import *

sys = make_system()
#for i in range(10, 20):
#    sys.gamma_list[i] = 1.0

print(sys.stability())

#gamma_lists = []
states = []

num_frames = 200
for i in range(num_frames):
    sys.step()
    #gamma_lists.append(sys.gamma_list.copy())
    states.append(sys.make_data())

ticks = sys.make_ticks()

# code based on:
# https://stackoverflow.com/questions/49165233/two-lines-matplotib-animation

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8))

ax1.set_xlim(0, sys.slice_length * sys.num_slices)
ax1.set_ylim(-0.1, 0.1)
line1, = ax1.plot([], [], color="r")

ax2.set_xlim(0, sys.slice_length * sys.num_slices)
ax2.set_ylim(-0.1, 0.1)
line2, = ax2.plot([], [], color="b")
line3, = ax2.plot([], [], color="y")

ax3.set_xlim(0, sys.slice_length * sys.num_slices)
ax3.set_ylim(-0.1, 0.1)
line4, = ax3.plot([], [], color="b")
line5, = ax3.plot([], [], color="y")


def update(n):
    (gamma, mu0_up, mu0_dn, mu0_hot_up, mu0_hot_dn) = states[n]
    line1.set_data(ticks, gamma)
    line2.set_data(ticks, mu0_up)
    line3.set_data(ticks, mu0_hot_up)
    line4.set_data(ticks, mu0_dn)
    line5.set_data(ticks, mu0_hot_dn)
    return [line1, line2, line3, line4, line5]


ani = animation.FuncAnimation(fig, update, num_frames, interval=100, blit=True)

plt.show()
