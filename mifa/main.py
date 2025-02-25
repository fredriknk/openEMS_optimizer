from mifa import ifa_simulation

def wifi():
    Sim_CSX = 'IFA.xml'
    showCad = True
    post_proc_only = True

    unit = 1e-3
    substrate_width = 21
    substrate_length = 30
    substrate_thickness = 1.5
    gndplane_position = 0
    substrate_cells = 4
    ifa_h = 6.965
    ifa_l = 21.019
    ifa_w1 = 2.961
    ifa_w2 = 1.235
    ifa_wf = 0.442
    ifa_fp = 4.086
    ifa_e = 0.5
    mifa_meander=1.5
    mifa_tipdistance=3.0
    mifa_meander_edge_distance=3.0
    substrate_epsR = 4.5
    feed_R = 50
    min_freq = 2.4e9
    center_freq = 2.45e9
    max_freq = 2.5e9
    min_size = 0.2 # minimum automesh size
    override_min_global_grid = 2 #none if not override
    max_timesteps = 800000
    plot = True

    #freq, s11_dB, Zin, P_in, hash_prefix = ifa_simulation(Sim_CSX=Sim_CSX,

def bt_groundplane():
    Sim_CSX = 'IFA.xml'
    showCad = True
    post_proc_only = False

    unit = 1e-3
    substrate_width = 21
    substrate_length = 83.1
    substrate_thickness = 1.5
    gndplane_position = -0.36
    substrate_cells = 4
    ifa_h = 6.072
    ifa_l = 17.9
    ifa_w1 = 0.613
    ifa_w2 = 0.4720
    ifa_wf = 0.425
    ifa_fp = 1.713
    ifa_e = 0.5
    mifa_meander=2
    mifa_tipdistance=3.0
    mifa_meander_edge_distance=3.0
    substrate_epsR = 4.5
    feed_R = 50
    min_freq = 2.4e9
    center_freq = 2.45e9
    max_freq = 2.5e9
    min_size = 0.201 # minimum automesh size
    override_min_global_grid = 2.4#None #none if not override
    max_timesteps = 100000
    plot = True

    freq, s11_dB, Zin, P_in, hash_prefix = ifa_simulation(Sim_CSX=Sim_CSX,
                                                    showCad=showCad,
                                                    post_proc_only=post_proc_only,
                                                    unit=post_proc_only,
                                                    substrate_width=substrate_width,
                                                    substrate_length=substrate_length,
                                                    substrate_thickness=substrate_thickness,
                                                    gndplane_position=gndplane_position,
                                                    substrate_cells=substrate_cells,
                                                    ifa_h=ifa_h,
                                                    ifa_l=ifa_l,
                                                    ifa_w1=ifa_w1,
                                                    ifa_w2=ifa_w2,
                                                    ifa_wf=ifa_wf,
                                                    ifa_fp=ifa_fp,
                                                    ifa_e=ifa_e,
                                                    mifa_meander=mifa_meander,
                                                    mifa_tipdistance=mifa_tipdistance,
                                                    mifa_meander_edge_distance=mifa_meander_edge_distance,
                                                    substrate_epsR=substrate_epsR,
                                                    feed_R=feed_R,
                                                    min_freq=min_freq,
                                                    center_freq=center_freq,
                                                    max_freq=max_freq,
                                                    min_size=min_size,
                                                    override_min_global_grid=override_min_global_grid,
                                                    max_timesteps=max_timesteps,
                                                    plot=plot)


def optimized_wifi():
    Sim_CSX = 'IFA.xml'
    showCad = True
    post_proc_only = False

    unit = 1e-3
    substrate_width = 21
    substrate_length = 40
    substrate_thickness = 1.5
    gndplane_position = -0.36
    substrate_cells = 4
    ifa_h = 5.5
    ifa_l = 21
    ifa_w1 = 1
    ifa_w2 = 1
    ifa_wf = 0.444
    ifa_fp = 4.069
    ifa_e = 0.5
    mifa_meander=2
    mifa_tipdistance=3.0
    mifa_meander_edge_distance=3.0
    substrate_epsR = 4.5
    feed_R = 50
    min_freq = 2.4e9
    center_freq = 2.45e9
    max_freq = 2.5e9
    min_size = 0.201 # minimum automesh size
    override_min_global_grid = 1.2 #none if not override
    max_timesteps = 100000
    plot = True

    freq, s11_dB, Zin, P_in, hash_prefix = ifa_simulation(Sim_CSX=Sim_CSX,
                                                    showCad=showCad,
                                                    post_proc_only=post_proc_only,
                                                    unit=post_proc_only,
                                                    substrate_width=substrate_width,
                                                    substrate_length=substrate_length,
                                                    substrate_thickness=substrate_thickness,
                                                    gndplane_position=gndplane_position,
                                                    substrate_cells=substrate_cells,
                                                    ifa_h=ifa_h,
                                                    ifa_l=ifa_l,
                                                    ifa_w1=ifa_w1,
                                                    ifa_w2=ifa_w2,
                                                    ifa_wf=ifa_wf,
                                                    ifa_fp=ifa_fp,
                                                    ifa_e=ifa_e,
                                                    mifa_meander=mifa_meander,
                                                    mifa_tipdistance=mifa_tipdistance,
                                                    mifa_meander_edge_distance=mifa_meander_edge_distance,
                                                    substrate_epsR=substrate_epsR,
                                                    feed_R=feed_R,
                                                    min_freq=min_freq,
                                                    center_freq=center_freq,
                                                    max_freq=max_freq,
                                                    min_size=min_size,
                                                    override_min_global_grid=override_min_global_grid,
                                                    max_timesteps=max_timesteps,
                                                    plot=plot)


def wifi_groundplane():
    Sim_CSX = 'IFA.xml'
    showCad = True
    post_proc_only = False

    unit = 1e-3
    substrate_width = 21
    substrate_length = 30
    substrate_thickness = 1.5
    gndplane_position = -0.36
    substrate_cells = 4
    ifa_h = 8.336
    ifa_l = 20.875
    ifa_w1 = 1.698
    ifa_w2 =  0.816
    ifa_wf = 1.398
    ifa_fp = 4.089
    ifa_e = 0.5
    mifa_meander=2
    mifa_tipdistance=2.0
    mifa_meander_edge_distance=2.0
    substrate_epsR = 4.5
    feed_R = 50
    min_freq = 2.4e9
    center_freq = 2.45e9
    max_freq = 2.5e9
    min_size = 0.2 # minimum automesh size
    override_min_global_grid = 1.2 #none if not override
    max_timesteps = 100000
    plot = True
    cleanup = True

    freq, s11_dB, Zin, P_in, hash_prefix = ifa_simulation(Sim_CSX=Sim_CSX,
                                                    showCad=showCad,
                                                    post_proc_only=post_proc_only,
                                                    unit=post_proc_only,
                                                    substrate_width=substrate_width,
                                                    substrate_length=substrate_length,
                                                    substrate_thickness=substrate_thickness,
                                                    gndplane_position=gndplane_position,
                                                    substrate_cells=substrate_cells,
                                                    ifa_h=ifa_h,
                                                    ifa_l=ifa_l,
                                                    ifa_w1=ifa_w1,
                                                    ifa_w2=ifa_w2,
                                                    ifa_wf=ifa_wf,
                                                    ifa_fp=ifa_fp,
                                                    ifa_e=ifa_e,
                                                    mifa_meander=mifa_meander,
                                                    mifa_tipdistance=mifa_tipdistance,
                                                    mifa_meander_edge_distance=mifa_meander_edge_distance,
                                                    substrate_epsR=substrate_epsR,
                                                    feed_R=feed_R,
                                                    min_freq=min_freq,
                                                    center_freq=center_freq,
                                                    max_freq=max_freq,
                                                    min_size=min_size,
                                                    override_min_global_grid=override_min_global_grid,
                                                    max_timesteps=max_timesteps,
                                                    plot=plot,
                                                    delete_simulation_files=cleanup)

def lte():
    #2024-12-22 18:43:45,044 - total seconds: 105.39, S11 at cf: -50.8712, Imp: 49.7R 0.1z, ifa_l: 140.389, ifa_h: 12.423, ifa_fp: 2.914, ifa_w1: 0.426, ifa_w2: 0.444, ifa_wf: 0.485, mifa_meander_edge_distance: 1.595, Res f: 0.830 GHz, S11 at res: -50.871, BW1: 0.79 GHz, BW2: 0.90 GHz, BW = 101.9 MHz - id 4235e8a4bce2
    Sim_CSX = 'IFA.xml'
    showCad = True
    post_proc_only = False
    unit = 1e-3
    substrate_width = 25.5
    substrate_length = 108
    substrate_thickness = 1.5
    gndplane_position = -0.36
    substrate_cells = 4
    ifa_h = 12.423
    ifa_l = 140.389
    ifa_w1 = 0.426
    ifa_w2 = 0.444
    ifa_wf = 0.485
    ifa_fp = 2.914
    ifa_e = 0.5
    mifa_meander=1.5
    md = 4
    mifa_tipdistance=1.5
    mifa_meander_edge_distance=1.595
    substrate_epsR = 4.5
    feed_R = 50
    min_freq = 0.791e9
    center_freq = 0.826e9
    max_freq = 0.862e9
    fc = 0.83e9-0.7e9
    min_size = 0.10 # minimum automesh size
    override_min_global_grid = 4 #none if not override
    max_timesteps = 2200000
    plot = True
    cleanup = False

    freq, s11_dB, Zin, P_in, hash_prefix = ifa_simulation(Sim_CSX=Sim_CSX,
                                                    showCad=showCad,
                                                    post_proc_only=post_proc_only,
                                                    unit=post_proc_only,
                                                    substrate_width=substrate_width,
                                                    substrate_length=substrate_length,
                                                    substrate_thickness=substrate_thickness,
                                                    gndplane_position=gndplane_position,
                                                    substrate_cells=substrate_cells,
                                                    ifa_h=ifa_h,
                                                    ifa_l=ifa_l,
                                                    ifa_w1=ifa_w1,
                                                    ifa_w2=ifa_w2,
                                                    ifa_wf=ifa_wf,
                                                    ifa_fp=ifa_fp,
                                                    ifa_e=ifa_e,
                                                    mifa_meander=mifa_meander,
                                                    mifa_tipdistance=mifa_tipdistance,
                                                    mifa_meander_edge_distance=mifa_meander_edge_distance,
                                                    substrate_epsR=substrate_epsR,
                                                    feed_R=feed_R,
                                                    fc=fc,
                                                    min_freq=min_freq,
                                                    center_freq=center_freq,
                                                    max_freq=max_freq,
                                                    min_size=min_size,
                                                    override_min_global_grid=override_min_global_grid,
                                                    max_timesteps=max_timesteps,
                                                    plot=plot)



#init main function
if __name__ == "__main__":
    lte()