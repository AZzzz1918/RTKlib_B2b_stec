#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot STEC and VTEC for G02, G05, G12 satellites from GAMP output.
Satellite ordering in output files: G01-G32, R01-R27, E01-E30, J01-J07, C01-C35
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os

# File paths
data_dir = r"E:\PPP\RTKLIB-B2b-master\RTKLIB-B2b-master_1\RTKLIB-B2b-master\1\result"
obs_name = "BDVL00AUS_R_20220150000_01D_30S_MO.22o"
stec_file = os.path.join(data_dir, obs_name + ".stec")
vtec_file = os.path.join(data_dir, obs_name + ".vtec")

# Satellite info: PRN -> 0-based index in the output array
# GPS: G01-G32 -> indices 0-31
SAT_INFO = {
    'G02': {'idx': 1},
    'G05': {'idx': 4},
    'G12': {'idx': 11},
}


def safe_float(s):
    """Parse float, handling -nan(ind) and other NaN notations from C printf output"""
    # C may output -nan(ind) which Python float() can't parse
    if 'nan' in s.lower():
        return np.nan
    try:
        val = float(s)
        if np.isnan(val) or np.isinf(val):
            return np.nan
        return val
    except (ValueError, OverflowError):
        return np.nan


def parse_stec(filepath, sat_info):
    """Parse STEC file, return dict of {sat_name: {'time': [], 'stec': []}}"""
    data = {name: {'time': [], 'stec': []} for name in sat_info}

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 8:
                continue

            # Parse time
            y, m, d, hh, mm, ss = [int(x) for x in fields[:6]]
            try:
                t = datetime(y, m, d, hh, mm, int(float(ss)))
            except ValueError:
                continue

            # Extract STEC for each requested satellite
            for name, info in sat_info.items():
                col = 8 + info['idx']
                if col < len(fields):
                    val = safe_float(fields[col])
                    # 99999.0 means no data
                    if abs(val - 99999.0) < 0.001:
                        val = np.nan
                    data[name]['time'].append(t)
                    data[name]['stec'].append(val)

    return data


def parse_vtec(filepath, sat_info):
    """Parse VTEC file, return dict of {sat_name: {'time': [], 'vtec': [], 'ipp_lat': [], 'ipp_lon': []}}"""
    data = {name: {'time': [], 'vtec': [], 'ipp_lat': [], 'ipp_lon': []} for name in sat_info}

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 8:
                continue

            # Parse time
            y, m, d, hh, mm, ss = [int(x) for x in fields[:6]]
            try:
                t = datetime(y, m, d, hh, mm, int(float(ss)))
            except ValueError:
                continue

            # Extract VTEC for each requested satellite
            # Format per satellite: vtec ipp_lon ipp_lat (3 columns)
            for name, info in sat_info.items():
                col = 8 + 3 * info['idx']
                if col + 2 < len(fields):
                    vtec = safe_float(fields[col])
                    ipp_lon = safe_float(fields[col + 1])
                    ipp_lat = safe_float(fields[col + 2])

                    # 99999.0 means no data
                    if abs(vtec - 99999.0) < 0.001:
                        vtec = np.nan
                    if abs(ipp_lon - 99999.0) < 0.001 or abs(ipp_lat - 99999.0) < 0.001:
                        ipp_lat = np.nan
                        ipp_lon = np.nan

                    data[name]['time'].append(t)
                    data[name]['vtec'].append(vtec)
                    data[name]['ipp_lat'].append(ipp_lat)
                    data[name]['ipp_lon'].append(ipp_lon)

    return data


def plot_stec(stec_data, sat_names):
    """Plot STEC time series for specified satellites"""
    fig, axes = plt.subplots(len(sat_names), 1, figsize=(14, 3.5 * len(sat_names)), sharex=True)
    if len(sat_names) == 1:
        axes = [axes]

    colors = ['#2196F3', '#FF5722', '#4CAF50']

    for ax, name, color in zip(axes, sat_names, colors):
        d = stec_data[name]
        times = d['time']
        stec = np.array(d['stec'])
        valid = ~np.isnan(stec)

        ax.plot(times, stec, '.', color=color, markersize=2, alpha=0.7, label=f'{name} STEC')
        ax.set_ylabel('STEC [TECU]', fontsize=12)
        ax.set_title(f'{name} — Slant TEC from GAMP PPP (BDVL station, 2022-01-15)', fontsize=13)
        ax.legend(loc='upper right', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.set_xlim(datetime(2022, 1, 15, 0, 0), datetime(2022, 1, 16, 0, 0))

    axes[-1].set_xlabel('Time [UTC]', fontsize=12)
    plt.tight_layout()
    return fig


def plot_vtec(vtec_data, sat_names):
    """Plot VTEC time series for specified satellites"""
    fig, axes = plt.subplots(len(sat_names), 1, figsize=(14, 3.5 * len(sat_names)), sharex=True)
    if len(sat_names) == 1:
        axes = [axes]

    colors = ['#2196F3', '#FF5722', '#4CAF50']

    for ax, name, color in zip(axes, sat_names, colors):
        d = vtec_data[name]
        times = d['time']
        vtec = np.array(d['vtec'])
        valid = ~np.isnan(vtec)

        ax.plot(times, vtec, '.', color=color, markersize=2, alpha=0.7, label=f'{name} VTEC')
        ax.set_ylabel('VTEC [TECU]', fontsize=12)
        ax.set_title(f'{name} — Vertical TEC from GAMP PPP (BDVL station, 2022-01-15)', fontsize=13)
        ax.legend(loc='upper right', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.set_xlim(datetime(2022, 1, 15, 0, 0), datetime(2022, 1, 16, 0, 0))

    axes[-1].set_xlabel('Time [UTC]', fontsize=12)
    plt.tight_layout()
    return fig


def plot_combined(stec_data, vtec_data, sat_names):
    """Combined STEC + VTEC plot, one row per satellite"""
    fig, axes = plt.subplots(len(sat_names), 2, figsize=(18, 3.5 * len(sat_names)))
    if len(sat_names) == 1:
        axes = axes.reshape(1, -1)

    colors = ['#2196F3', '#FF5722', '#4CAF50']

    for i, (name, color) in enumerate(zip(sat_names, colors)):
        # STEC subplot
        ax_stec = axes[i, 0]
        d_stec = stec_data[name]
        stec = np.array(d_stec['stec'])
        ax_stec.plot(d_stec['time'], stec, '.', color=color, markersize=2, alpha=0.7)
        ax_stec.set_ylabel('STEC [TECU]', fontsize=12)
        ax_stec.set_title(f'{name} STEC', fontsize=13)
        ax_stec.grid(True, alpha=0.3)
        ax_stec.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax_stec.set_xlim(datetime(2022, 1, 15, 0, 0), datetime(2022, 1, 16, 0, 0))

        # VTEC subplot
        ax_vtec = axes[i, 1]
        d_vtec = vtec_data[name]
        vtec = np.array(d_vtec['vtec'])
        ax_vtec.plot(d_vtec['time'], vtec, '.', color=color, markersize=2, alpha=0.7)
        ax_vtec.set_ylabel('VTEC [TECU]', fontsize=12)
        ax_vtec.set_title(f'{name} VTEC', fontsize=13)
        ax_vtec.grid(True, alpha=0.3)
        ax_vtec.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax_vtec.set_xlim(datetime(2022, 1, 15, 0, 0), datetime(2022, 1, 16, 0, 0))

    for ax in axes[-1]:
        ax.set_xlabel('Time [UTC]', fontsize=12)

    fig.suptitle('GAMP PPP: STEC & VTEC — BDVL station, 2022-01-15 (DOY 015)', fontsize=15, fontweight='bold')
    plt.tight_layout()
    return fig


def main():
    sat_names = ['G02', 'G05', 'G12']

    print(f"Parsing STEC file: {stec_file}")
    stec_data = parse_stec(stec_file, SAT_INFO)
    for name in sat_names:
        n_valid = np.sum(~np.isnan(stec_data[name]['stec']))
        mn = np.nanmin(stec_data[name]['stec']) if n_valid > 0 else float('nan')
        mx = np.nanmax(stec_data[name]['stec']) if n_valid > 0 else float('nan')
        print(f"  {name}: {len(stec_data[name]['time'])} epochs, {n_valid} valid STEC, "
              f"range=[{mn:.2f}, {mx:.2f}] TECU")

    print(f"Parsing VTEC file: {vtec_file}")
    vtec_data = parse_vtec(vtec_file, SAT_INFO)
    for name in sat_names:
        n_valid = np.sum(~np.isnan(vtec_data[name]['vtec']))
        mn = np.nanmin(vtec_data[name]['vtec']) if n_valid > 0 else float('nan')
        mx = np.nanmax(vtec_data[name]['vtec']) if n_valid > 0 else float('nan')
        print(f"  {name}: {len(vtec_data[name]['time'])} epochs, {n_valid} valid VTEC, "
              f"range=[{mn:.2f}, {mx:.2f}] TECU")

    # Generate combined plot (STEC + VTEC side by side)
    print("Generating combined STEC+VTEC plot...")
    fig_comb = plot_combined(stec_data, vtec_data, sat_names)
    outpath_comb = os.path.join(data_dir, "G02_G05_G12_STEC_VTEC.png")
    fig_comb.savefig(outpath_comb, dpi=150, bbox_inches='tight')
    print(f"Saved: {outpath_comb}")

    # Generate STEC-only plot
    print("Generating STEC plot...")
    fig_stec = plot_stec(stec_data, sat_names)
    outpath_stec = os.path.join(data_dir, "G02_G05_G12_STEC.png")
    fig_stec.savefig(outpath_stec, dpi=150, bbox_inches='tight')
    print(f"Saved: {outpath_stec}")

    # Generate VTEC-only plot
    print("Generating VTEC plot...")
    fig_vtec = plot_vtec(vtec_data, sat_names)
    outpath_vtec = os.path.join(data_dir, "G02_G05_G12_VTEC.png")
    fig_vtec.savefig(outpath_vtec, dpi=150, bbox_inches='tight')
    print(f"Saved: {outpath_vtec}")

    print("All plots generated successfully!")


if __name__ == '__main__':
    main()
