# openEMStesting

Testing the OpenEMS FDTD solver for antenna optimization.

The entire simulation stack is opensource and free to use for anyone commercial or private, allowing you to create matched antennas for your designs.  

The Mifa_groundplane.py is an automated IFA/MIFA generator which allows for automated optimization of antenna design

The antenna generation is optimized for the avaliable groundplane so in essence you can just design a pcb, input the groundplane parameters and the script will generate a ifa/mifa for you based on the avalible space on the pcb.

Avaliable paramaters are as follows:

```plaintext 
#############################################################
#|------------- substrate_width -------------------|
# _______________________________________________      _ substrate_thickness
#| A  ifa_e      |----------ifa_l(total length)-| |\   \-gndplane_position 
#| V____          _______________     __________| | |   \_0 point
#|               |    ___  ___   |___|  ______  | | |
#|         ifa_h |   |   ||   |_________|    |  | | |_ mifa_meander_edge_distance 
#|               |   |   ||  mifa_meander    |__| | |_ mifa_tipdistance
#|               |   |   ||                   w2  | | |                  
#|_______________|___|___||_______________________| |_|
#| <---ifa_e---->| w1|   wf\                      | |
#|               |__fp___|  \                     | |
#|                       |    feed point          | |
#|                       |                        | | substrate_length
#|<- substrate_width/2 ->|                        | |
#|                                                | |
#|________________________________________________| |
# \________________________________________________\|
#############################################################
# Note: It's not checked whether your settings make sense, so check
#       graphical output carefully.
```

The antenna generator makes an IFA antenna if ifa_l < substrate_length, if it is longer it will meander down to mifa_tipdistance, and if it is even longer than it will start to create meanders with a minimum distance of mifa_meander_edge_distance going back towards the feedpoint.

there arent any checks and balances, so please check the "extreme cases" if you use the optimization script

The antenna generator also meshes the antenna elements and groundplane using the 2/3 1/3 rule, so each element will get a centerline and a box of size min_size going 1/3 the element and 1/3 outside it. A

If gndplane_position < 0, then the script will move the groundplane down the specified distance and place a via at the short ciruit stub

## Example

The parameters used

```Python
Sim_CSX = 'IFA.xml'
showCad = True
post_proc_only = True

unit = 1e-3
substrate_width = 21
substrate_length = 20
substrate_thickness = 1.5
gndplane_position = -0.36
substrate_cells = 4
ifa_h = 5.500
ifa_l = 40#IM changing this variable
ifa_w1 = 1
ifa_w2 = 1
ifa_wf = 1
ifa_fp = 5.000
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
max_timesteps = 800000
plot = True
```

Im changing the ifa_l variable keeping everything else the same.

<table>
  <tr>
    <td align="center"><strong>15mm < substrate_width </strong><br>
      <img src="pictures/SHORT.PNG" alt="Short" width="200"/>
    </td>
    <td align="center"><strong>23mm > substrate_width</strong><br>
      <img src="pictures/MEDIUM.PNG" alt="Medium" width="200"/>
    </td>
    <td align="center"><strong>40mm >> substrate_width</strong><br>
      <img src="pictures/LONG.PNG" alt="Long" width="200"/>
    </td>
  </tr>
</table>

## Optimizer

The optimizer is a differential evolution multivariable optimizer, its much faster with single variable gradient optimizers, but i get much better results with multivariable optimizers though they do take a lot longer to finish.

### 2.4 GHZ
Here is an example of the optimizer output after a 6 hour simulation for a 21x30mm substrate

<table>
  <tr>
    <td align="center"><strong>S11 chart </strong><br>
      <img src="pictures/2_4GHz_optimized.png" alt="Short" height="300px"/>
    </td>
    <td align="center"><strong>Antenna</strong><br>
      <img src="pictures/Optimized.PNG" alt="Medium" height="300px"/>
    </td>
    </tr>
</table>


## Installation

You need to follow the instructions on the [openEMS website](https://docs.openems.de/install/package.html#install-readymade-windows-package-src) to install the software, and then you can run the script by running the following command in the terminal:

I would reccomend unzipping to this location, then you install the CSXCAD and openEMS python packages using the following commands:

```bash
cd C:\opt\openEMS\python
pip install CSXCAD-0.6.2-cp310-cp310-win_amd64.whl
pip install openEMS-0.0.33-cp310-cp310-win_amd64.whl
```

You may need to make an .env file in the same directory as the script with the following parameters:

```bash
# openEMS ENV file
# ----------------
rootdir = "C:/opt/openEMS"
csxcad_location = "C:/opt/openEMS/AppCSXCAD"
```
