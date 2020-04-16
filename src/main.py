#import matplotlib.pyplot as plt
#plt.plot([1, 2, 3, 4], [1, 8, 27, 64])
#plt.ylabel('some numbers')
#plt.show()

from src.system import *

sys = make_system()
sys.gamma_list[40] = 1.0

print(sys.stability())

for i in range(10):
    for j in range(50):
        sys.step()
    sys.plot()