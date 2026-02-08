import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 14})

sizes = np.array([
    4,
    16,
    64,
    256,
    4096,
    16384,
    65536,
    262144,
    1048576
]) * 4 # message size in bytes

times = np.array([
    1e-6,
    2e-6,
    2e-6,
    2e-6,
    7e-6,
    1.4e-5,
    5.7e-5,
    1.04e-4,
    3.73e-4
])

def comm_time(m, alpha, beta):
    return alpha + m * beta

filtered_sizes = sizes[2:]
filtered_times = times[2:]

breakpoint = 4

smaller_sizes = filtered_sizes[:breakpoint]
smaller_times = filtered_times[:breakpoint]

bigger_sizes = filtered_sizes[breakpoint:]
bigger_times = filtered_times[breakpoint:]

popt1, pcov1 = curve_fit(comm_time, smaller_sizes, smaller_times)
popt2, pcov2 = curve_fit(comm_time, bigger_sizes, bigger_times)

plt.semilogy(filtered_sizes, filtered_times, 'o')
plt.semilogy(smaller_sizes, comm_time(smaller_sizes, *popt1), '-', label=rf'Fit 1: {popt1[0]:.2e}$\cdot m +${popt1[1]:.2e}')
plt.semilogy(bigger_sizes, comm_time(bigger_sizes, *popt2), '-', label=rf'Fit 2: {popt2[0]:.2e}$\cdot m +${popt2[1]:.2e}')

plt.xlabel('Message size (bytes)')
plt.ylabel('Time (s)')
plt.title("Two process Ping-Pong (different nodes)")
plt.legend()
plt.tight_layout()
plt.grid(alpha=.2)
plt.savefig('./pingpong_double_node.pdf')