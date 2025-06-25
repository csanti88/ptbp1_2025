import os
import shutil
import subprocess
#use Rmats filtering .py file
#adapted from https://github.com/Xinglab/rmats-turbo-tutorial/tree/main
# https://github.com/Xinglab/rmats-turbo-tutorial/blob/main/scripts/rmats_filtering.py

# Path to rmats_filtering.py
sc_rmats_filter = r"..\rmats_filtering.py"

# rMATS results path
rmats_input_dir = r"H:\bw\20250527_wtE16\new2_join_wt28_vs_e16"

# Output directory
output_dir = r"H:\bw\20250527_wtE16\new2_join_wt28_vs_e16_tmp_output"
os.makedirs(output_dir, exist_ok=True)

# Event and count types
event_array = ["SE", "MXE", "RI", "A5SS", "A3SS"]
counttype_array = ["JC", "JCEC"]

# Filtering parameters

# readCov = 20
# minPSI = 0.05
# maxPSI = 0.95
# sigFDR = 0.01
# sigDeltaPSI = 0.05
# bgFDR = 0.5
# bgWithinGroupDeltaPSI = 0.2

params = "0,0.01,0.99,0.05,0.10,0.5,0.2"

for event in event_array:
    for counttype in counttype_array:
        input_file = os.path.join(rmats_input_dir, f"{event}.MATS.{counttype}.txt")
        if os.path.exists(input_file):
            print(f"Copying: {input_file}")
            shutil.copy(input_file, output_dir)
            copied_file = os.path.join(output_dir, f"{event}.MATS.{counttype}.txt")
            print(f"Running filter on: {copied_file}")
            subprocess.run([
                "python", sc_rmats_filter,
                copied_file,
                params
            ])
        else:
            print(f"File not found: {input_file}")