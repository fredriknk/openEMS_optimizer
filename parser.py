import pandas as pd
import re

# Define the log file path
log_file_path = 'c:/Users/fnk/pythonprojects/openEMS/logs/diffevolution_log.txt'

# Read the log file
with open(log_file_path, 'r') as file:
    log_data = file.read()

# Define the regex pattern to extract the data
pattern = re.compile(
    r'(?P<timestamp>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}) - total seconds: (?P<total_seconds>[\d.]+), '
    r'ifa_l: (?P<ifa_l>[\d.]+), ifa_h: (?P<ifa_h>[\d.]+), ifa_fp: (?P<ifa_fp>[\d.]+), '
    r'ifa_w1: (?P<ifa_w1>[\d.]+), ifa_w2: (?P<ifa_w2>[\d.]+), ifa_wf: (?P<ifa_wf>[\d.]+), '
    r'S11 at cf: (?P<s11_cf>[-\d.]+), Imp: (?P<imp>[^,]+), Res f: (?P<res_f>[\d.]+) GHz, '
    r'S11 at res: (?P<s11_res>[-\d.]+), BW1: (?P<bw1>[\d.]+) GHz, BW2: (?P<bw2>[\d.]+) GHz, '
    r'BW = (?P<bw>[\d.]+) MHz - id (?P<id>\w+)'
)

# Find all matches in the log data
matches = pattern.findall(log_data)

# Convert matches to a DataFrame
df = pd.DataFrame(matches, columns=[
    'timestamp', 'total_seconds', 'ifa_l', 'ifa_h', 'ifa_fp', 'ifa_w1', 'ifa_w2', 'ifa_wf',
    's11_cf', 'imp', 'res_f', 's11_res', 'bw1', 'bw2', 'bw', 'id'
])

# Convert appropriate columns to numeric types
numeric_columns = [
    'total_seconds', 'ifa_l', 'ifa_h', 'ifa_fp', 'ifa_w1', 'ifa_w2', 'ifa_wf',
    's11_cf', 'res_f', 's11_res', 'bw1', 'bw2', 'bw'
]
for col in numeric_columns:
    df[col] = pd.to_numeric(df[col])

# Sort by 's11_cf' and select the top 2 entries with the lowest values
best_antennas = df.sort_values(by='s11_cf').head(10)

# Display the best antennas
print(best_antennas)