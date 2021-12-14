import matplotlib.pyplot as plt
import kshell_utilities as ksutil

def nodes_4():
    res_128 = ksutil.get_timing_data("4_nodes/128")
    res_256 = ksutil.get_timing_data("4_nodes/256")

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

def nodes_1():
    res_128 = ksutil.get_timing_data("1_node/128")
    res_256 = ksutil.get_timing_data("1_node/256")

    fig, ax = plt.subplots()
    ax.bar(["128 cores", "256 virtual cores (SMT)"], [res_128, res_256])
    ax.set_title("V50 GXPF1A\n1 node, 256 virtual cores per node")
    ax.set_ylabel("Time [s]")
    plt.show()

    fig, ax = plt.subplots()
    ax.bar(["128 cores", "256 virtual cores (SMT)"], [res_128/res_128, res_256/res_128])
    ax.set_title("V50 GXPF1A\n1 nodes, 256 virtual cores per node")
    ax.set_ylabel("Relative time")
    plt.show()


if __name__ == "__main__":
    # nodes_4()
    nodes_1()