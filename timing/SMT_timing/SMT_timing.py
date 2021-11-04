import matplotlib.pyplot as plt
import kshell_utilities as ksutil

def main():
    res_128 = ksutil.get_timing_data("128")
    res_256 = ksutil.get_timing_data("256")

    fig, ax = plt.subplots()
    ax.bar(["128 cores", "256 virtual cores (SMT)"], [res_128, res_256])
    ax.set_title("V50 GXPF1A\n4 nodes, 256 virtual cores per node")
    ax.set_ylabel("Time [s]")
    plt.show()

    fig, ax = plt.subplots()
    ax.bar(["128 cores", "256 virtual cores (SMT)"], [res_128/res_128, res_256/res_128])
    ax.set_title("V50 GXPF1A\n4 nodes, 256 virtual cores per node")
    ax.set_ylabel("Relative time")
    plt.show()


if __name__ == "__main__":
    main()