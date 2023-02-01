GS_FREE_PROTON = 5.585
GS_FREE_NEUTRON = -3.826

recommended_quenching_factors = {
    "GCLSTsdpfsdgix5pn.snt": f"0.75*GS_FREE = {round(0.75*GS_FREE_PROTON, 3), round(0.75*GS_FREE_NEUTRON, 3)}",
    "gs8.snt": f"0.75*GS_FREE = {round(0.75*GS_FREE_PROTON, 3), round(0.75*GS_FREE_NEUTRON, 3)}",
    "jun45.snt": f"0.7*GS_FREE = {round(0.7*GS_FREE_PROTON, 3), round(0.7*GS_FREE_NEUTRON, 3)}",
    "gxpf1a.snt": f"0.9*GS_FREE = {round(0.9*GS_FREE_PROTON, 3), round(0.9*GS_FREE_NEUTRON, 3)}",
    "gxpf1.snt": f"0.9*GS_FREE = {round(0.9*GS_FREE_PROTON, 3), round(0.9*GS_FREE_NEUTRON, 3)}",
    "sdpf-mu.snt": f"0.9*GS_FREE = {round(0.9*GS_FREE_PROTON, 3), round(0.9*GS_FREE_NEUTRON, 3)}",
    "sn100pn.snt": f"0.7*GS_FREE = {round(0.7*GS_FREE_PROTON, 3), round(0.7*GS_FREE_NEUTRON, 3)}"
}
