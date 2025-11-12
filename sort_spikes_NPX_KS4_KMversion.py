# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 10:02:10 2024

Based on Isabel's sort spikes, but updated a bit for the changes in spikeinterface

I have not run this on the cluster yet, so I might split into preprocessing, kilosort, and postprocessing steps


"""

import spikeinterface.full as si
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from spikeinterface.sortingcomponents.peak_localization import localize_peaks
from spikeinterface.exporters import export_to_phy
from probeinterface import ProbeGroup, write_prb

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from glob import glob
import os
import time
import pandas as pd
import fastparquet

from kilosort import run_kilosort
from kilosort import io

import argparse

## Set paths 
# lets run some more on gpu
#base_folder = Path('/n/data2/hms/neurobio/harvey/Kevin/Videos/KM33')
#spikeglx_folder = base_folder / '250216_g0'

#base_folder = Path('/n/data2/hms/neurobio/harvey/Kevin/Videos/KM33-34')
#spikeglx_folder = base_folder / '250218_g0'

#base_folder = Path('/n/data2/hms/neurobio/harvey/Kevin/Videos/KM37')
#spikeglx_folder = base_folder / '250216_g0'

base_folder = Path('/n/data2/hms/neurobio/harvey/Kevin/Videos/KM37-38')
spikeglx_folder = base_folder / '250219_g0'

## Arg parser to get from arguments
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mouseName",type=str)
parser.add_argument("-d", "--recordingDate",type=str)
args = parser.parse_args()

#base_folder = Path('/n/data2/hms/neurobio/harvey/Kevin/Videos/' + args.mouseName)
base_folder = Path('/n/scratch/users/b/biw842/' + args.mouseName)
spikeglx_folder = base_folder / args.recordingDate

#data_folder = Path('/n/scratch/users/k/km315/' + args.mouseName)
#spikeglx_folder_load = data_folder / args.recordingDate


## Set control

do_preprocess  = 1
do_kilosort    = 1
do_postprocess = 1
save_figs      = 0
save_preprocess = 0

## Run stuff
# folders to save things
figure_path = os.path.join(spikeglx_folder, 'figures')
if not os.path.exists(figure_path): os.makedirs(figure_path)
output_path = os.path.join(spikeglx_folder, 'Spike_Sorting')
if not os.path.exists(output_path): os.makedirs(output_path)
probe_path = os.path.join(spikeglx_folder, 'preprocess_spikeinterface','prb.prb') # what's this?
if not os.path.exists(os.path.join(spikeglx_folder, 'preprocess_spikeinterface')): os.makedirs(os.path.join(spikeglx_folder, 'preprocess_spikeinterface'))

print('Sorting spikes for: ')
print(spikeglx_folder)

## Load recording stream
stream_names, stream_ids = si.get_neo_streams('spikeglx', spikeglx_folder) # would be load if doing separate locations
raw_rec = si.read_spikeglx(spikeglx_folder, stream_name = 'imec0.ap', load_sync_channel=False)


## TEMP ONLY ONCE FOR KM33-34 split recording!
#print('Concatenating with - ')
#spikeglx_folder_2 = base_folder / '250427-2_g0'
#raw_rec_2 = si.read_spikeglx(spikeglx_folder_2, stream_name = 'imec0.ap', load_sync_channel=False)

#raw_rec = si.concatenate_recordings([raw_rec, raw_rec_2])
#print(spikeglx_folder_2)

#print('Concatenating very specific recordings should not see this if not expecting this!!')

#output_path = os.path.join(spikeglx_folder, 'Spike_Sorting_concatenated')
#if not os.path.exists(output_path): os.makedirs(output_path)

## PREPROCESSING STEPS
# Honestly this is pretty quick so I prob should just do it every time...



if do_preprocess:
    # a) high pass filter (300hz? 400 hz?)
    rec_hp = si.highpass_filter(raw_rec, freq_min=300.)
    # b) phase shift <- what does this one do?
    rec_shift = si.phase_shift(rec_hp)
    # c) remove bad channels (how does it know?)
    bad_channel_ids, channel_labels = si.detect_bad_channels(rec_shift)
    rec_rm_bad = si.interpolate_bad_channels(rec_shift, bad_channel_ids)
    print('bad_channel_ids', bad_channel_ids)
    # d) common average referecing <- ok this one i know
    rec_cmr = si.common_reference(rec_rm_bad, operator='median', reference='global')

else:
    try: 
        rec_cmr = si.BinaryFolderRecording.load_from_folder(os.path.join(spikeglx_folder,'preprocess_spikeinterface'))
        print('Loading preprocessed data from binary')
    except:
        print('you need to preprocess!')
        
if save_figs:
    print('Plotting different preprocess steps')
    fig, axs = plt.subplots(ncols=3, figsize=(15, 10))
    si.plot_timeseries(rec_hp, backend='matplotlib',  clim=(-50, 50), ax=axs[0])
    si.plot_timeseries(rec_shift, backend='matplotlib',  clim=(-50, 50), ax=axs[1])
    si.plot_timeseries(rec_cmr, backend='matplotlib',  clim=(-50, 50), ax=axs[2])
    for i, label in enumerate(('filter', 'phase_shift', 'cmr')):
        axs[i].set_title(label)
    plt.savefig(os.path.join(figure_path,'preprocessing.png'))
    
    #Plot subset of traces 
    print('Plotting subset of traces')
    plt_kwargs = {'linewidth':0.05}
    fig, ax = plt.subplots(1,1, figsize=(10,20),**plt_kwargs)
    channels = rec_cmr.get_channel_ids()[::20]
    si.plot_timeseries(rec_cmr, time_range=[200,205],figure=fig,channel_ids=channels,backend='matplotlib',mode='line')
    ax.set_yticks(np.arange(1,len(channels)+1))
    ax.set_yticklabels(np.arange(1,len(channels)+1))
    ax.set_xticklabels(np.arange(6))
    ax.set_ylabel('Channels');ax.set_xlabel('Time(s)')
    plt.axis('off')
    plt.savefig(os.path.join(figure_path,'raw_filtered_timeseries_abbrev.pdf'),format='pdf')
    
    #Plot full set of traces 
    print('Plotting full traces')
    plt_kwargs = {'linewidth':0.05}
    fig, ax = plt.subplots(1,1, figsize=(10,200),**plt_kwargs)
    channels = rec_cmr.get_channel_ids()
    si.plot_timeseries(rec_cmr, time_range=[11,16],figure=fig,channel_ids=channels,backend='matplotlib',mode='line')
    ax.set_yticks(np.arange(1,len(channels)+1))
    #label channel index 
    ax.set_yticklabels(np.arange(1,len(channels)+1))
    ax.set_xticklabels(np.arange(6))
    ax.set_ylabel('Channels');ax.set_xlabel('Time(s)')
    plt.axis('off')
    plt.savefig(os.path.join(figure_path,'raw_filtered_timeseries.pdf'),format='pdf')
        
    # Plot noise levels
    print('Plotting noise')
    noise_levels_microV = si.get_noise_levels(rec_cmr, return_scaled=True)
    noise_levels_int16 = si.get_noise_levels(rec_cmr, return_scaled=False)
    
    fig, ax = plt.subplots()
    _ = ax.hist(noise_levels_microV, bins=np.arange(5, 30, 2.5))
    ax.set_xlabel('noise  [microV]')
    ax.set_ylabel('occurance (a.u.)')
    plt.axis('off')
    plt.savefig(os.path.join(figure_path, 'noise_levels.pdf'), format='pdf')
    plt.savefig(os.path.join(figure_path, 'noise_levels.png'))
    
    # Plot drift
    print('Plotting drift')
    job_kwargs = dict(n_jobs=40, chunk_duration='1s', progress_bar=True)
    peaks = detect_peaks(rec_cmr,  method='locally_exclusive', noise_levels=noise_levels_int16,
                     detect_threshold=5, radius_um=50., **job_kwargs)
    
    peak_locations = localize_peaks(rec_cmr, peaks, method='center_of_mass', radius_um=50., **job_kwargs)
    fs = rec_cmr.sampling_frequency
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(peaks['sample_index'] / fs, peak_locations['y'], color='k', marker='.',  alpha=0.002)
    plt.axis('off')
    plt.savefig(os.path.join(figure_path, 'drift.pdf'), format='pdf',dpi=50)
    plt.savefig(os.path.join(figure_path, 'drift.png'))
    
    # Plot probe map
    fig, ax = plt.subplots(figsize=(15,10))
    probe_ylower = 700
    probe_yupper = 2200
    si.plot_probe_map(rec_cmr, ax=axs[2], with_channel_ids=False)
    axs[2].set_ylim(probe_ylower, probe_yupper)
    axs[2].scatter(peak_locations['x'], peak_locations['y'], color='purple', alpha=0.002)
    plt.axis('off')
    plt.savefig(os.path.join(figure_path, 'probe_map.pdf'), format='pdf',dpi=50)
    plt.savefig(os.path.join(figure_path, 'probe_map.png'))
    
if save_preprocess:
    #Save preprocessed data to binary file 
    job_kwargs = dict(n_jobs=40, chunk_duration='1s', progress_bar=True,chunk_memory="10M")
    print('Saving to binary')
    rec_cmr = rec_cmr.save(folder=os.path.join(spikeglx_folder,'preprocess_spikeinterface'),overwrite=True, format='binary', **job_kwargs)
    print('Saved!')
    
## KILOSORT


if do_kilosort:
    # First load the probe data. This is better than using a default?
    probe = raw_rec.get_probe()
    pg = ProbeGroup()
    pg.add_probe(probe)
    write_prb(probe_path, pg)
    probe = io.load_probe(probe_path)
    fs = raw_rec.get_sampling_frequency()
    
    # Set paths for spike sorting
    settings = {'fs': fs, 'n_chan_bin': 384,'data_dir':output_path } # should this be 385?
    #rec_cmr = si.BinaryFolderRecording.load_from_folder(os.path.join(spikeglx_folder,'preprocess_spikeinterface'))
    print(rec_cmr) # this shoudl be geenrated from above?
    wrapper = io.RecordingExtractorAsArray(rec_cmr)
    sorting_folder = output_path + '/kilosort4'
    
    try:
        sorting_KS4 = si.read_kilosort(folder_path=sorting_folder)
        print('Loading pre-computed sorting')
    except:
        start = time.time()
        print ('Launching kilosort4...')
        ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = run_kilosort(
        settings=settings, probe=probe, filename=spikeglx_folder, file_object=wrapper,
            results_dir=sorting_folder)
        print('sorting finished in ',time.time() - start )
        sorting_KS4 = si.read_kilosort(folder_path=sorting_folder)

else:
    probe = io.load_probe(probe_path)
    sorting_folder = output_path + '/kilosort4'
    sorting_KS4 = si.read_kilosort(folder_path=sorting_folder)
    
# POST PROCESSING
# this will be a bit different from isabels since ill use the sorting analyzer
if do_postprocess:
    #save the primary channels for each identified unit 
    templates = np.load(os.path.join(sorting_folder,'templates.npy'))
    channel_map = np.load(os.path.join(sorting_folder,'channel_map.npy'))
    peak_chan_idx = np.squeeze(np.argmax(np.max(templates,1) - np.min(templates,1),1))
    channel_map = np.squeeze(channel_map)
    peak_channels = np.squeeze(channel_map[peak_chan_idx])
    np.save(os.path.join(output_path,'primary_channels.npy'),peak_channels)
    
    # extract waveforms:
    job_kwargs = dict(n_jobs=40, chunk_duration='1s', progress_bar=True) # this can't be 40 on local machien
    waveforms_path = sorting_folder + '/waveforms_ks4'
    
    # create postprocessing analyzer
    # **??** was I supposed to put job_kwargs here?
    sorting_analyzer = si.create_sorting_analyzer(sorting=sorting_KS4, recording=rec_cmr, folder=waveforms_path, sparse=True, **job_kwargs)
    # compute some important information about units
    sorting_analyzer.compute('random_spikes', method='uniform', max_spikes_per_unit=500,**job_kwargs)
    sorting_analyzer.compute('waveforms', ms_before=1.5, ms_after=2.5, **job_kwargs)
    sorting_analyzer.compute('principal_components', n_components=3, mode='by_channel_local')
    sorting_analyzer.compute('templates', operators = ['average','median','std'])
    sorting_analyzer.compute('noise_levels')
    
    # save quality metrics
    qc = sorting_analyzer.compute("quality_metrics").get_data() # IDK this might be the replacement
    
    qc = pd.DataFrame(qc)
    qc_path = output_path + '/quality_metrics.parquet'
    qc.to_parquet(qc_path, compression='gzip')
    
    #curate using quality metrics 
    #good_units = qc[(qc.snr>4) & (qc.isi_violations_ratio<1)&(qc.firing_rate>0.01)]
    good_units = qc[(qc.snr>2) & (qc.isi_violations_ratio<1)&(qc.firing_rate>0.05)&(qc.nn_hit_rate>0.5)&(qc.amplitude_cutoff<0.1)&(qc.presence_ratio>0.9)]
    good_units = np.array(good_units.index)
    np.savez(qc_path, good_units=good_units)
    print('Saved units and quality control metrics...')
    print('You have', len(good_units), 'good units!')

    # save sorting analyzer?
    # - useful if want to plot waveforms, etc
    print('Saving analyzer')
    good_sorting_analyzer = sorting_analyzer.select_units(unit_ids=sorting_analyzer.unit_ids[good_units], format="binary_folder", folder=waveforms_path)
    good_sorting_analyzer.save_as(folder=os.path.join(waveforms_path, 'good_sorting_analyzer.zarr'), format='zarr')

    # export to phy (slow)
    phy_path = os.path.join(output_path,'phy')
    if os.path.exists(phy_path): 
        print('')
    print('Exporting to phy')
    export_to_phy(sorting_analyzer, output_path + '/phy', **job_kwargs)
    
    print('Sorting complete!')
