import matplotlib.pyplot as plt
import matplotlib.animation as animation

from src.system import *

sys = make_system()
#for i in range(10, 20):
#    sys.gamma_list[i] = 1.0

print(sys.stability())

gamma_lists = []

num_frames = 200
for i in range(num_frames):
    sys.step()
    gamma_lists.append(sys.gamma_list.copy())

ticks = sys.make_ticks()

# code based on:
# https://stackoverflow.com/questions/49165233/two-lines-matplotib-animation

fig, ax = plt.subplots()
ax.set_xlim(0, sys.slice_length * sys.num_slices)
ax.set_ylim(0, 1)
line1, = ax.plot([], [], color="r")


def update(n):
    line1.set_data(ticks, gamma_lists[n])
    return [line1]


ani = animation.FuncAnimation(fig, update, num_frames, interval=100, blit=True)

plt.show()
