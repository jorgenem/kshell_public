import os
import numpy as np
import matplotlib.pyplot as plt
import kshell_utilities as ksutil

def get_filenames(name_beginning: str):
    res = [i for i in os.listdir() if i.startswith(name_beginning)]
    res.sort(key=lambda fname: int(fname.split("=")[-1]))
    
    return res

def main():
    """
    Plot timing and memory comparisons for different combinations of
    n_block and n_nodes.
    """
    si28_fnames = get_filenames("Si28")
    v47_fnames = get_filenames("V47")
    v50_4nodes_fnames = get_filenames("V50_n_nodes=04")
    v50_8nodes_fnames = get_filenames("V50_n_nodes=08")
    v50_16nodes_fnames = get_filenames("V50_n_nodes=16")

    si28_m_dim = ksutil.count_dim(
        model_space_filename = f"{si28_fnames[0]}/usda.snt",
        partition_filename = f"{si28_fnames[0]}/Si28_usda_p.ptn",
        print_dimensions = False
    )[1][-1]

    v47_m_dim = ksutil.count_dim(
        model_space_filename = f"{v47_fnames[0]}/gxpf1a.snt",
        partition_filename = f"{v47_fnames[0]}/V47_gxpf1a_n.ptn",
        print_dimensions = False
    )[1][-1]
    
    v50_4nodes_m_dim = ksutil.count_dim(
        model_space_filename = f"{v50_4nodes_fnames[0]}/gxpf1a.snt",
        partition_filename = f"{v50_4nodes_fnames[0]}/V50_gxpf1a_p.ptn",
        print_dimensions = False
    )[1][-1]

    v50_8nodes_m_dim = ksutil.count_dim(
        model_space_filename = f"{v50_8nodes_fnames[0]}/gxpf1a.snt",
        partition_filename = f"{v50_8nodes_fnames[0]}/V50_gxpf1a_p.ptn",
        print_dimensions = False
    )[1][-1]

    v50_16nodes_m_dim = ksutil.count_dim(
        model_space_filename = f"{v50_16nodes_fnames[0]}/gxpf1a.snt",
        partition_filename = f"{v50_16nodes_fnames[0]}/V50_gxpf1a_p.ptn",
        print_dimensions = False
    )[1][-1]

    n_si28_fnames = len(si28_fnames)
    n_v47_fnames = len(v47_fnames)
    n_v50_4nodes_fnames = len(v50_4nodes_fnames)
    n_v50_8nodes_fnames = len(v50_8nodes_fnames)
    n_v50_16nodes_fnames = len(v50_16nodes_fnames)
    
    si28_times = np.zeros(n_si28_fnames)
    si28_memory = np.zeros(n_si28_fnames)
    si28_tr_times = np.zeros(n_si28_fnames)
    si28_n_block = np.zeros(n_si28_fnames)
    
    v47_times = np.zeros(n_v47_fnames)
    v47_memory = np.zeros(n_v47_fnames)
    v47_tr_times = np.zeros(n_v47_fnames)
    v47_n_block = np.zeros(n_v47_fnames)
    
    v50_4nodes_times = np.zeros(n_v50_4nodes_fnames)
    v50_4nodes_memory = np.zeros(n_v50_4nodes_fnames)
    v50_4nodes_tr_times = np.zeros(n_v50_4nodes_fnames)
    v50_4nodes_n_block = np.zeros(n_v50_4nodes_fnames)
    
    v50_8nodes_times = np.zeros(n_v50_8nodes_fnames)
    v50_8nodes_memory = np.zeros(n_v50_8nodes_fnames)
    v50_8nodes_tr_times = np.zeros(n_v50_8nodes_fnames)
    v50_8nodes_n_block = np.zeros(n_v50_8nodes_fnames)
    
    v50_16nodes_times = np.zeros(n_v50_16nodes_fnames)
    v50_16nodes_memory = np.zeros(n_v50_16nodes_fnames)
    v50_16nodes_tr_times = np.zeros(n_v50_16nodes_fnames)
    v50_16nodes_n_block = np.zeros(n_v50_16nodes_fnames)

    for i in range(n_si28_fnames):
        si28_times[i] = ksutil.get_timing_data(f"{si28_fnames[i]}/log_Si28_usda_m0p.txt")
        si28_memory[i] = ksutil.get_memory_usage(f"{si28_fnames[i]}/log_Si28_usda_m0p.txt")
        si28_tr_times[i] = ksutil.get_timing_data(f"{si28_fnames[i]}/log_Si28_usda_tr_m0p_m0p.txt")
        si28_n_block[i] = int(si28_fnames[i].split("=")[-1])

    for i in range(n_v47_fnames):
        v47_times[i] = ksutil.get_timing_data(f"{v47_fnames[i]}/log_V47_gxpf1a_m1n.txt")
        v47_memory[i] = ksutil.get_memory_usage(f"{v47_fnames[i]}/log_V47_gxpf1a_m1n.txt")
        v47_tr_times[i] = ksutil.get_timing_data(f"{v47_fnames[i]}/log_V47_gxpf1a_tr_m1n_m1n.txt")
        v47_n_block[i] = int(v47_fnames[i].split("=")[-1])

    for i in range(n_v50_4nodes_fnames):
        v50_4nodes_times[i] = ksutil.get_timing_data(f"{v50_4nodes_fnames[i]}/log_V50_gxpf1a_m0p.txt")
        v50_4nodes_memory[i] = ksutil.get_memory_usage(f"{v50_4nodes_fnames[i]}/log_V50_gxpf1a_m0p.txt")
        v50_4nodes_tr_times[i] = ksutil.get_timing_data(f"{v50_4nodes_fnames[i]}/log_V50_gxpf1a_tr_m0p_m0p.txt")
        v50_4nodes_n_block[i] = int(v50_4nodes_fnames[i].split("=")[-1])

    for i in range(n_v50_8nodes_fnames):
        v50_8nodes_times[i] = ksutil.get_timing_data(f"{v50_8nodes_fnames[i]}/log_V50_gxpf1a_m0p.txt")
        v50_8nodes_memory[i] = ksutil.get_memory_usage(f"{v50_8nodes_fnames[i]}/log_V50_gxpf1a_m0p.txt")
        v50_8nodes_tr_times[i] = ksutil.get_timing_data(f"{v50_8nodes_fnames[i]}/log_V50_gxpf1a_tr_m0p_m0p.txt")
        v50_8nodes_n_block[i] = int(v50_8nodes_fnames[i].split("=")[-1])

    for i in range(n_v50_16nodes_fnames):
        v50_16nodes_times[i] = ksutil.get_timing_data(f"{v50_16nodes_fnames[i]}/log_V50_gxpf1a_m0p.txt")
        v50_16nodes_memory[i] = ksutil.get_memory_usage(f"{v50_16nodes_fnames[i]}/log_V50_gxpf1a_m0p.txt")
        v50_16nodes_tr_times[i] = ksutil.get_timing_data(f"{v50_16nodes_fnames[i]}/log_V50_gxpf1a_tr_m0p_m0p.txt")
        v50_16nodes_n_block[i] = int(v50_16nodes_fnames[i].split("=")[-1])

    fig, ax = plt.subplots()
    ax.plot([4, 8, 16], np.array([v50_4nodes_times[0], v50_8nodes_times[0], v50_16nodes_times[0]])/v50_4nodes_times[0], "--.", color="darkred")
    ax.set_title(r"kshell only, $^{50}$V" + "\nn_block=0")
    ax.set_xlabel("n_nodes")
    ax.set_ylabel("Relative time")
    plt.show()

    si28_times /= si28_times[0]
    si28_tr_times /= si28_tr_times[0]
    v47_times /= v47_times[0]
    v47_tr_times /= v47_tr_times[0]
    v50_4nodes_times /= v50_4nodes_times[0]
    v50_4nodes_tr_times /= v50_4nodes_tr_times[0]
    v50_8nodes_times /= v50_8nodes_times[0]
    v50_8nodes_tr_times /= v50_8nodes_tr_times[0]
    v50_16nodes_times /= v50_16nodes_times[0]
    v50_16nodes_tr_times /= v50_16nodes_tr_times[0]

    fig, ax = plt.subplots()
    ax.plot(si28_n_block, si28_times, "--.", color="black", label=f"M-dim = {si28_m_dim}, n_nodes = 1")
    ax.plot(v47_n_block, v47_times, "--.", color="grey", label=f"M-dim = {v47_m_dim}, n_nodes = 1")
    ax.plot(v50_4nodes_n_block, v50_4nodes_times, "--.", color="darkred", label=f"M-dim = {v50_4nodes_m_dim}, n_nodes = 4")
    ax.plot(v50_8nodes_n_block, v50_8nodes_times, "--.", color="darkblue", label=f"M-dim = {v50_8nodes_m_dim}, n_nodes = 8")
    ax.plot(v50_16nodes_n_block, v50_16nodes_times, "--.", color="darkgreen", label=f"M-dim = {v50_16nodes_m_dim}, n_nodes = 16")
    ax.set_title(r"kshell only, $^{50}$V")
    ax.set_xlabel("n_block")
    ax.set_ylabel("Relative time")
    ax.set_xticks(range(0, 20+1, 2))
    ax.legend()
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(si28_n_block, si28_tr_times, "--.", color="black", label=f"M-dim = {si28_m_dim}, n_nodes = 1")
    ax.plot(v47_n_block, v47_tr_times, "--.", color="grey", label=f"M-dim = {v47_m_dim}, n_nodes = 1")
    ax.plot(v50_4nodes_n_block, v50_4nodes_tr_times, "--.", color="darkred", label=f"M-dim = {v50_4nodes_m_dim}, n_nodes = 4")
    ax.plot(v50_8nodes_n_block, v50_8nodes_tr_times, "--.", color="darkblue", label=f"M-dim = {v50_8nodes_m_dim}, n_nodes = 8")
    ax.plot(v50_16nodes_n_block, v50_16nodes_tr_times, "--.", color="darkgreen", label=f"M-dim = {v50_16nodes_m_dim}, n_nodes = 16")
    ax.set_title(r"transit only, $^{50}$V")
    ax.set_xlabel("n_block")
    ax.set_ylabel("Relative time")
    ax.set_xticks(range(0, 20+1, 2))
    ax.legend()
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(si28_n_block, si28_memory, "--.", color="black", label=f"M-dim = {si28_m_dim}, n_nodes = 1")
    ax.plot(v47_n_block, v47_memory, "--.", color="grey", label=f"M-dim = {v47_m_dim}, n_nodes = 1")
    ax.plot(v50_4nodes_n_block, v50_4nodes_memory/4, "--.", color="darkred", label=f"M-dim = {v50_4nodes_m_dim}, n_nodes = 4")
    ax.plot(v50_8nodes_n_block, v50_8nodes_memory/8, "--.", color="darkblue", label=f"M-dim = {v50_8nodes_m_dim}, n_nodes = 8")
    ax.plot(v50_16nodes_n_block, v50_16nodes_memory/16, "--.", color="darkgreen", label=f"M-dim = {v50_16nodes_m_dim}, n_nodes = 16")
    ax.set_title("kshell")
    ax.set_xlabel("n_block")
    ax.set_ylabel("Memory usage per node [GB]")
    ax.set_xticks(range(0, 20+1, 2))
    ax.legend()
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(v50_16nodes_n_block, v50_16nodes_memory/v50_16nodes_memory[0], "--.", color="darkgreen", label=f"M-dim = {v50_16nodes_m_dim}, n_nodes = 16")
    ax.set_title("kshell")
    ax.set_xlabel("n_block")
    ax.set_ylabel("Normalized memory usage")
    ax.set_xticks(range(0, 20+1, 2))
    ax.legend()
    plt.show()

if __name__ == "__main__":
    main()