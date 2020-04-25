import matplotlib.pyplot as plt
import matplotlib.animation as animation

from src.system import *

sys = make_system()
states = []

num_frames = 200
for i in range(num_frames):
    print(i / num_frames)
    for j in range(1):
        sys.step()
    states.append(sys.make_data())

(ticks_slices, ticks_planes) = sys.make_ticks()

# code based on:
# https://stackoverflow.com/questions/49165233/two-lines-matplotib-animation

fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, figsize=(10, 8))

ax1.set_xlim(0, sys.slice_length * sys.num_slices)
ax1.set_ylim(-0.2, 0.1)
ax1.set_ylabel("gamma\n(eV)")
ax1.axhline(y=0, color="g")
line1, = ax1.plot([], [], color="r")

ax2.set_xlim(0, sys.slice_length * sys.num_slices)
ax2.set_ylim(-0.1, 0.1)
ax2.set_ylabel("mu_0 up\n(eV)")
ax2.axhline(y=0, color='g')
line2, = ax2.plot([], [], color="b")
#line3, = ax2.plot([], [], color="y")

ax3.set_xlim(0, sys.slice_length * sys.num_slices)
ax3.set_ylim(-0.1, 0.1)
ax3.set_ylabel("mu_0 dn\n(eV)")
ax3.axhline(y=0, color='g')
line4, = ax3.plot([], [], color="b")
#line5, = ax3.plot([], [], color="y")

ax4.set_xlim(0, sys.slice_length * sys.num_slices)
ax4.set_ylabel("N_hot\n(nm^-3)")
ax4.set_ylim(0.0, 3.0)
line6, = ax4.plot([], [], color="b")
line7, = ax4.plot([], [], color="y")

ax5.set_xlim(0, sys.slice_length * sys.num_slices)
ax5.set_ylim(-1.0, 1.0)
ax5.set_ylabel("J_hot\n(nm^-2 fs^-1)")
ax5.axhline(y=0, color='g')
line8, = ax5.plot([], [], color="b")
line9, = ax5.plot([], [], color="y")

ax6.set_xlim(0, sys.slice_length * sys.num_slices)
ax6.set_ylim(-0.5, 1.0)
ax6.set_ylabel("J_spin\n(nm^2 fs^-1)")
ax6.axhline(y=0, color='g')
line10, = ax6.plot([], [], color="r")


def init():
    line1.set_data([], [])
    line2.set_data([], [])
    #line3.set_data([], [])
    line4.set_data([], [])
    #line5.set_data([], [])
    line6.set_data([], [])
    line7.set_data([], [])
    line8.set_data([], [])
    line9.set_data([], [])
    line10.set_data([], [])
    #return [line1, line2, line3, line4, line5, line6, line7, line8, line9, line10]
    return [line1, line2, line4, line6, line7, line8, line9, line10]


def update(n):
    (gamma, mu0_up, mu0_dn, hot_up, hot_dn, mu0_hot_up, mu0_hot_dn, j_hot_up, j_hot_dn,
     j_up, j_dn, j_spin) = states[n]
    line1.set_data(ticks_slices, gamma)
    line2.set_data(ticks_slices, mu0_up)
    #line3.set_data(ticks_slices, mu0_hot_up)
    line4.set_data(ticks_slices, mu0_dn)
    #line5.set_data(ticks_slices, mu0_hot_dn)
    line6.set_data(ticks_slices, hot_up)
    line7.set_data(ticks_slices, hot_dn)
    line8.set_data(ticks_planes, j_hot_up)
    line9.set_data(ticks_planes, j_hot_dn)
    line10.set_data(ticks_planes, j_spin)
    #return [line1, line2, line3, line4, line5, line6, line7, line8, line9, line10]
    return [line1, line2, line4, line6, line7, line8, line9, line10]


anim = animation.FuncAnimation(fig, update, init_func=init, frames=num_frames, interval=100, blit=True)

show = True
if show:
    plt.show()
else:
    writer = animation.FFMpegWriter(fps=30.0)
    anim.save("output/animation.mp4", writer=writer)
