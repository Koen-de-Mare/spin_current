import matplotlib.pyplot as plt
import matplotlib.animation as animation

from src.system import *

sys = make_system()
states = []

num_frames = 200
for i in range(num_frames):
    print(i / num_frames)
    for j in range(2):
        sys.step()
    states.append(sys.make_data())

(ticks_slices, ticks_planes) = sys.make_ticks()

# code based on:
# https://stackoverflow.com/questions/49165233/two-lines-matplotib-animation

fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(10, 8))

ax1.set_xlim(0, sys.slice_length * sys.num_slices)
ax1.set_ylim(-0.1, 0.1)
ax1.axhline(y=0, color='g')
line1, = ax1.plot([], [], color="r")

ax2.set_xlim(0, sys.slice_length * sys.num_slices)
ax2.set_ylim(-0.1, 0.1)
ax2.axhline(y=0, color='g')
line2, = ax2.plot([], [], color="b")
line3, = ax2.plot([], [], color="y")

ax3.set_xlim(0, sys.slice_length * sys.num_slices)
ax3.set_ylim(-0.1, 0.1)
ax3.axhline(y=0, color='g')
line4, = ax3.plot([], [], color="b")
line5, = ax3.plot([], [], color="y")

ax4.set_xlim(0, sys.slice_length * sys.num_slices)
ax4.set_ylim(0.0, 10.0)

line6, = ax4.plot([], [], color="b")
line7, = ax4.plot([], [], color="y")

ax5.set_xlim(0, sys.slice_length * sys.num_slices)
ax5.set_ylim(-1.0, 1.0)
ax5.axhline(y=0, color='g')
line8, = ax5.plot([], [], color="r")

def update(n):
    (gamma, mu0_up, mu0_dn, hot_up, hot_dn, mu0_hot_up, mu0_hot_dn, j_hot, j_up, j_dn) = states[n]
    line1.set_data(ticks_slices, gamma)
    line2.set_data(ticks_slices, mu0_up)
    line3.set_data(ticks_slices, mu0_hot_up)
    line4.set_data(ticks_slices, mu0_dn)
    line5.set_data(ticks_slices, mu0_hot_dn)
    line6.set_data(ticks_slices, hot_up)
    line7.set_data(ticks_slices, hot_dn)
    line8.set_data(ticks_planes, j_hot)
    return [line1, line2, line3, line4, line5, line6, line7, line8]

anim = animation.FuncAnimation(fig, update, num_frames, interval=100, blit=True)

show = True
if show:
    plt.show()
else:
    writer = animation.FFMpegWriter(fps=30.0)
    anim.save("output/animation.mp4", writer=writer)
