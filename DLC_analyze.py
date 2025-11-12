# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 10:02:10 2024

Based on Isabel's sort spikes, but updated a bit for the changes in spikeinterface

I have not run this on the cluster yet, so I might split into preprocessing, kilosort, and postprocessing steps


"""

import deeplabcut as dlc
import argparse



## Arg parser to get from arguments
parser = argparse.ArgumentParser()

parser.add_argument("-d", "--recordingDate",type=str)
args = parser.parse_args()

#base_folder = Path('/n/data2/hms/neurobio/harvey/Kevin/Videos/' + args.mouseName)
base_folder = Path('/n/scratch/users/b/biw842/' + args.recordingDate + '.mp4')



cfg = r"/n/scratch/users/b/biw842/config.yaml"
vid = base_folder
print(args.recordingDate)
# Re-run with auto_track so DLC uses the tracker settings you edited in inference_cfg.yaml
dlc.analyze_videos(
    cfg, vid,
    engine=dlc.Engine.PYTORCH,
    batch_size=8,        # raise if VRAM allows (e.g., 32)
    auto_track=True,
    n_tracks=2,           # set to your real # of animals
    save_as_csv=True,
)
print("finished")