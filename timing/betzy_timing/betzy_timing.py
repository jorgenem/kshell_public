import matplotlib.pyplot as plt

def betzy_timing():
    V50_10_states_n_nodes = [4, 8, 16, 32, 64]
    V50_10_states_timing = [32.27, 22.58, 27.129, 50.83, 16.19]
    plt.title(r"$^{50}$V 10 states")
    plt.plot(V50_10_states_n_nodes, V50_10_states_timing, "--.")
    plt.xlabel("Nodes")
    plt.ylabel("Time [s]")
    plt.show()

    V50_400_states_timing = [861.95107, 509.98084, 329.2985, 269.06916, 199.01395]
    V50_400_states_n_nodes = [4, 8, 16, 32, 64]
    plt.title(r"$^{50}$V 400 states")
    plt.plot(V50_400_states_n_nodes, V50_400_states_timing, "--.")
    plt.xlabel("Nodes")
    plt.ylabel("Time [s]")
    plt.show()

    V50_400_states_timing = [156.28706, 135.44382, 70.85086]
    V50_400_states_M_scheme_dim = [2003318, 1112120, 134526]
    plt.title(r"$^{50}$V 400 states")
    plt.plot(V50_400_states_M_scheme_dim, V50_400_states_timing, "--.")
    plt.xlabel("M-scheme dim for spin 0")
    plt.ylabel("Time [s]")
    plt.show()
    
if __name__ == "__main__":
    betzy_timing()