import json
from GaAntenna import ga_simulation
# JSON string (replace this with your actual JSON string)
json_data = """{"total_seconds":61.75,"fitness":"-14.77012","S11_at_cf":["-12.726","-17.142"],"params":{"Sim_CSX":"IFA.xml","unit":0.001,"substrate_width":25.5,"substrate_length":107,"substrate_thickness":1.5,"substrate_epsR":4.5,"gndplane_position":0,"substrate_cells":4,"ant_h":45,"ant_l":24.5,"ant_fp":5,"ant_e":0.5,"feed_R":50,"min_freq":800000000.0,"center_freq":1300000000.0,"max_freq":1800000000.0,"overlap":0.1,"fc":700000000.0,"max_timesteps":10000,"override_min_global_grid":null,"plot":false,"showCad":false,"post_proc_only":false,"delete_simulation_files":true,"antenna_grid":[[0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,1],[1,1,1,1,0,0,0,1,1,0,1,0,0,0,0,0,0,1,1,1],[0,1,1,0,1,1,1,0,1,0,0,1,1,0,0,0,0,1,0,1],[0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,1,0,1],[0,0,0,1,1,1,1,1,0,1,0,0,0,1,0,0,0,1,1,0],[1,1,1,1,0,1,0,0,0,1,0,0,0,0,0,1,0,0,1,1],[0,1,1,1,0,1,0,1,0,0,0,1,1,1,0,0,0,1,1,1],[1,0,0,0,1,0,1,1,0,1,0,1,0,0,1,0,1,1,0,1],[0,1,1,1,1,1,1,0,0,0,0,0,0,1,0,1,1,1,1,0],[0,0,1,1,1,0,1,0,1,1,1,0,1,0,1,1,0,0,0,1],[0,1,0,0,1,1,1,1,0,1,1,0,0,0,1,0,1,1,0,0],[0,1,0,1,0,1,0,1,1,0,1,1,1,0,1,1,0,1,1,0],[0,1,0,1,0,0,1,0,0,0,1,0,0,1,1,1,0,0,1,1],[0,0,1,0,1,1,1,1,1,0,1,1,1,0,0,0,0,0,1,1],[0,1,1,0,1,1,1,1,1,0,1,0,1,1,0,1,0,1,1,0],[1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,0,1,0,1],[0,1,0,0,0,0,1,1,0,0,1,0,0,0,1,1,1,0,1,0],[1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,0,0,0,0],[1,1,0,1,1,1,0,1,1,0,1,1,1,0,0,1,0,0,0,1],[1,1,1,1,0,0,1,0,0,0,1,1,0,1,0,0,1,0,0,1],[1,1,1,0,1,1,0,1,1,0,0,1,0,0,0,0,1,1,1,1],[1,0,0,0,0,0,1,1,0,1,0,0,1,1,0,0,0,1,1,0],[1,1,1,0,1,0,0,1,1,0,0,0,1,1,0,1,1,1,0,0],[1,1,0,0,1,1,1,0,0,0,1,1,1,0,1,1,1,0,1,1],[0,1,1,1,1,1,0,0,1,0,0,0,0,1,0,0,1,1,0,1],[0,1,0,1,1,0,1,0,0,0,1,1,0,1,0,0,1,0,0,1],[1,1,0,0,0,1,1,1,0,1,0,0,0,1,0,1,1,0,0,0],[1,1,1,1,0,0,1,1,1,0,1,1,1,0,1,0,1,0,1,0],[1,0,0,1,0,1,1,1,0,1,0,1,0,1,1,0,0,1,1,1],[0,1,1,0,0,1,1,1,0,1,1,1,0,1,0,0,0,1,0,1],[0,0,1,1,1,0,1,0,0,0,1,1,0,0,0,1,0,0,1,0],[0,1,1,1,1,0,1,1,1,1,0,1,0,1,1,1,0,1,1,1],[1,0,1,1,1,0,0,1,1,0,0,1,1,1,1,1,0,0,0,1],[0,1,0,0,1,0,1,1,0,0,1,1,0,1,1,1,1,0,1,1],[0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,1,0],[1,0,1,0,1,1,0,1,1,0,1,1,0,0,0,1,0,1,0,1],[0,1,0,0,1,1,1,1,1,0,0,0,1,0,1,1,0,0,0,1],[0,1,1,0,1,0,1,1,1,0,0,1,1,0,1,1,0,0,0,1],[1,0,0,1,1,0,1,0,1,0,1,1,1,1,0,0,1,1,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],"randomhash":496813,"frequencies":[830000000.0,1800000000.0],"numThreads":5,"lambdamultiplier":2,"cell_size_x":1.225,"cell_size_y":1.125,"feed_point":[-6.7375,53.5625,1.5]}}
"""
# Parse the JSON data
parsed_data = json.loads(json_data)

# Access individual components
total_seconds = parsed_data["total_seconds"]
s11_at_cf = parsed_data["S11_at_cf"]
params = parsed_data["params"]

# Example of accessing the `antenna_grid` as a list of lists
antenna_grid = params["antenna_grid"]
print("Antenna grid (first row):", antenna_grid[0])

# Print other parameters or data as needed
print("Total seconds:", total_seconds)
print("S11 at center frequency:", s11_at_cf)
print("Feed point:", params["feed_point"])
import numpy as np

def deserialize_data(data):
    if isinstance(data, list):
        # Try to convert lists back to np.ndarray, if possible
        try:
            return np.array(data)  # Convert list to numpy array if possible
        except:
            return [deserialize_data(item) for item in data]  # Handle nested lists recursively
    elif isinstance(data, dict):
        # Recursively handle dictionaries
        return {key: deserialize_data(value) for key, value in data.items()}
    else:
        # Base case: return the item itself if it's not a list or dict
        return data


params["plot"] = True  # Set the plot parameter to True
params["showCad"] = True  # Set the showCad parameter to True
params["post_proc_only"] = False  # Set the post_proc_only parameter to False
params["antenna_grid"] = np.array(deserialize_data(antenna_grid))  # Update the antenna grid parameter
params["max_timesteps"] = 120000  # Update the max_timesteps parameter
#params["numthreads"] = 0  # Add a new parameter numthreads
#params['mesh_divide']= 3 # Add a new parameter mesh_divide
#params["override_min_global_grid"] = 0.25  # Set the override_min_global_grid parameter to None
ga_simulation(parameters=params)  # Call the simulation function with the parsed parameters