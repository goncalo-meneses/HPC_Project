import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 14})

def comm_time(m, alpha, beta):
    return alpha + m * beta

sizes = np.array([2 ** i for i in range(2, 21, 2)]) * 4    # message size in bytes
times = np.array([
    0.000000,
    0.000012,
    0.000004,
    0.000005,
    0.000006,
    0.000011,
    0.000035,
    0.000075,
    0.000156,
    0.000615
])

filtered_sizes = sizes[2:]
filtered_times = times[2:]

breakpoint = 4

smaller_sizes = filtered_sizes[:breakpoint]
smaller_times = filtered_times[:breakpoint]

bigger_sizes = filtered_sizes[breakpoint:]
bigger_times = filtered_times[breakpoint:]

popt1, pcov1 = curve_fit(comm_time, smaller_sizes, smaller_times)
popt2, pcov2 = curve_fit(comm_time, bigger_sizes, bigger_times)

plt.plot(filtered_sizes, filtered_times, 'o')
plt.plot(smaller_sizes, comm_time(smaller_sizes, *popt1), '-', label=rf'Fit 1: {popt1[0]:.2e}$\cdot m +${popt1[1]:.2e}')
plt.plot(bigger_sizes, comm_time(bigger_sizes, *popt2), '-', label=rf'Fit 2: {popt2[0]:.2e}$\cdot m +${popt2[1]:.2e}')

plt.xlabel('Message size (bytes)')
plt.ylabel('Time (s)')
plt.title("Two process Ping-Pong (same node)")
plt.legend()
plt.tight_layout()
plt.grid(alpha=.2)
plt.savefig('./pingpong_single_node.pdf')