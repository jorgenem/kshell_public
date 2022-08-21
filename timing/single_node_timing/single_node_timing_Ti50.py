import numpy as np
import matplotlib.pyplot as plt
import kshell_utilities as ksutil

def main():
    """
    Produce timing plots for single node performance using data from
    Ti50.
    """
    paths = ["Ti50_no_trunc", "Ti50_trunc_1", "Ti50_trunc_2"]
    N = len(paths)

    mdim_arr = np.zeros(N)
    kshell_time = np.zeros(N)
    transit_time = np.zeros(N)

    for i in range(N):
        kshell_time[i] = \
            ksutil.get_timing_data(path=f"{paths[i]}/log_Ti50_gxpf1a_m0p.txt")
        transit_time[i] = \
            ksutil.get_timing_data(path=f"{paths[i]}/log_Ti50_gxpf1a_tr_m0p_m0p.txt")
        M, mdim, jdim = ksutil.count_dim(
            model_space_filename = f"{paths[i]}/gxpf1a.snt",
            partition_filename = f"{paths[i]}/Ti50_gxpf1a_p.ptn",
            print_dimensions = False,
            debug = False
        )
        mdim_arr[i] = mdim[-1]

    total_time = kshell_time + transit_time

    kshell_slope = (mdim_arr[2] - mdim_arr[1])/(kshell_time[2] - kshell_time[1])
    transit_slope = (mdim_arr[2] - mdim_arr[1])/(transit_time[2] - transit_time[1])
    total_slope = (mdim_arr[2] - mdim_arr[1])/(total_time[2] - total_time[1])

    plt.plot(kshell_time, mdim_arr, ".--", color="black", label=f"KSHELL, slope: {kshell_slope:.0f}")
    plt.plot(transit_time, mdim_arr, ".--", color="darkgrey", label=f"transit: {transit_slope:.0f}")
    plt.plot(total_time, mdim_arr, ".--", color="grey", label=f"total: {total_slope:.0f}")
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel("M-scheme dimension for lowest spin")
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()