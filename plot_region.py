#!/usr/bin/env python3
"""
Genomic Region Methylation Visualization

This script creates genome-browser style plots showing methylation levels 
at non-B DNA motifs within specific genomic regions.

Author: Your Name
License: MIT
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from matplotlib.patches import Rectangle

def plot_genomic_region_methylation(bed_files, chromosome, start, end, output_dir="./plots"):
    """
    Plot methylation levels in a specific genomic region - Genome Browser style
    
    Parameters:
    bed_files (dict or list): Dictionary with labels as keys and file paths as values,
                             or list of file paths
    chromosome (str): Chromosome to analyze (e.g., 'chrX_MATERNAL')
    start (int): Start position
    end (int): End position
    output_dir (str): Directory to save output plots
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # BED file column names for BEDMethyl format
    columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 
               'thickEnd', 'itemRgb', 'coverage', 'methylation_pct', 'methylated', 
               'unmethylated', 'context1', 'context2', 'context3', 'context4', 'context5']
    
    # Handle both dict and list input formats
    if isinstance(bed_files, dict):
        files_to_process = bed_files
    elif isinstance(bed_files, list):
        # Create labels from filenames
        files_to_process = {os.path.basename(f).replace('.bed', ''): f for f in bed_files}
    else:
        raise ValueError("bed_files must be a dictionary or list")
    
    # Filter to only include existing files
    existing_files = {label: path for label, path in files_to_process.items() 
                     if path and os.path.exists(path)}
    
    if not existing_files:
        print("No valid files found!")
        return
    
    # Create the plot
    fig, axes = plt.subplots(len(existing_files), 1, figsize=(15, 3*len(existing_files)), sharex=True)
    if len(existing_files) == 1:
        axes = [axes]
    
    # Define colors for different modification types
    mod_colors = {
        '5mC': '#1f77b4',      # blue
        '5hmC': '#d62728',     # red
        'G4': '#9467bd',       # purple
        'GQ': '#9467bd',       # purple
        'Z': '#ff7f0e',        # orange
        'STR': '#2ca02c',      # green
        'IR': '#e377c2',       # pink
        'MR': '#8c564b',       # brown
        'DR': '#7f7f7f',       # gray
        'default': '#17becf'   # cyan
    }
    
    # Process each file
    for i, (label, file_path) in enumerate(existing_files.items()):
        print(f"Processing {label}...")
        
        try:
            # Load data
            df = pd.read_csv(file_path, sep='\t', header=None, names=columns)
            
            # Filter for the specific chromosome and region
            region_data = df[(df['chrom'] == chromosome) & 
                           (df['start'] >= start) & (df['end'] <= end)]
            
            print(f"  Found {len(region_data)} sites in region {chromosome}:{start}-{end}")
            
            if len(region_data) == 0:
                axes[i].text(0.5, 0.5, f'No data for {label} in region', 
                           transform=axes[i].transAxes, ha='center', va='center')
                axes[i].set_ylabel(label, rotation=0, labelpad=50, ha='right', va='center')
                continue
            
            # Extract clean motif and modification names for display
            display_name = extract_clean_name(label)
            
            # Determine color based on label
            color = mod_colors['default']
            label_lower = label.lower()
            
            # Priority: modification type first, then motif type
            for mod_key, mod_color in mod_colors.items():
                if mod_key.lower() in label_lower and mod_key != 'default':
                    color = mod_color
                    break
            
            # Plot each methylation site as a rectangle
            for _, row in region_data.iterrows():
                # Rectangle parameters
                rect_start = row['start']
                rect_width = max(100, row['end'] - row['start'])  # Minimum width of 100bp
                rect_height = row['methylation_pct']
                
                # Create rectangle with methylation height and color based on type
                rect = Rectangle((rect_start, 0), rect_width, rect_height,
                               facecolor=color, alpha=0.7, edgecolor='black', linewidth=0.3)
                axes[i].add_patch(rect)
            
            # Customize the subplot with clean display name
            axes[i].set_ylim(0, 100)
            axes[i].set_ylabel(f'{display_name}\nMethylation %', rotation=0, labelpad=50, ha='right', va='center')
            axes[i].grid(True, alpha=0.3, axis='y')
            axes[i].set_xlim(start, end)
            
            # Add sample size info
            sample_size = len(region_data)
            axes[i].text(0.02, 0.95, f'n={sample_size}', transform=axes[i].transAxes,
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                        fontsize=8, ha='left', va='top')
            
        except Exception as e:
            print(f"Error processing {label}: {e}")
            display_name = extract_clean_name(label)
            axes[i].text(0.5, 0.5, f'Error loading {display_name}', 
                        transform=axes[i].transAxes, ha='center', va='center')
            axes[i].set_ylabel(display_name, rotation=0, labelpad=50, ha='right', va='center')
    
    # Set x-axis labels
    axes[-1].set_xlabel(f'Genomic Position ({chromosome})')
    axes[-1].set_xlim(start, end)
    
    # Format x-axis to show genomic positions nicely
    axes[-1].ticklabel_format(style='plain', axis='x')
    
    # Set title
    plt.suptitle(f'Methylation Levels in {chromosome}:{start:,}-{end:,}', fontsize=14, y=0.98)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # Save plot
    region_str = f"{chromosome}_{start}_{end}".replace(':', '_')
    output_file = os.path.join(output_dir, f'genomic_region_methylation_{region_str}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved to: {output_file}")

def extract_clean_name(label):
    """
    Extract clean display name from file label
    e.g., 'GQ_5mC_intersect_file' -> 'G-Quadruplex 5mC'
    """
    # Motif mapping for clean names
    motif_mapping = {
        'G4': 'G-Quadruplex',
        'GQ': 'G-Quadruplex',
        'IR': 'Inverted Repeat',
        'MR': 'Mirror Repeat',
        'DR': 'Direct Repeat',
        'STR': 'Short Tandem Repeat',
        'Z': 'Z-DNA',
        'SLIPPED': 'Slipped Structure',
        'CRUCIFORM': 'Cruciform',
        'TRIPLEX': 'Triplex'
    }
    
    label_upper = label.upper()
    
    # Extract modification type
    mod_type = ""
    if "5HMC" in label_upper:
        mod_type = "5hmC"
    elif "5MC" in label_upper:
        mod_type = "5mC"
    else:
        # Look for any methylation pattern
        if "MC" in label_upper:
            mod_type = "mC"
    
    # Extract motif type
    motif_type = "Unknown"
    for motif_key, display_name in motif_mapping.items():
        if motif_key.upper() in label_upper:
            motif_type = display_name
            break
    
    # If still unknown, try simple extraction
    if motif_type == "Unknown":
        # Try to get motif from label
        parts = label.split('_')
        for part in parts:
            if part.upper() in motif_mapping:
                motif_type = motif_mapping[part.upper()]
                break
        if motif_type == "Unknown" and len(parts) > 0:
            motif_type = parts[0].title()
    
    # Combine motif and modification type
    if mod_type:
        return f"{motif_type} {mod_type}"
    else:
        return motif_type

def main():
    parser = argparse.ArgumentParser(description='Plot methylation levels in a genomic region')
    parser.add_argument('--bed-files', nargs='+', required=True,
                       help='Input BED files with methylation data')
    parser.add_argument('--chromosome', required=True,
                       help='Chromosome to analyze (e.g., chrX_MATERNAL)')
    parser.add_argument('--start', type=int, required=True,
                       help='Start position')
    parser.add_argument('--end', type=int, required=True,
                       help='End position')
    parser.add_argument('--output-dir', default='./plots',
                       help='Directory to save output plots')
    
    args = parser.parse_args()
    
    # Create bed_files dictionary from command line arguments
    bed_files = args.bed_files
    
    plot_genomic_region_methylation(
        bed_files=bed_files,
        chromosome=args.chromosome,
        start=args.start,
        end=args.end,
        output_dir=args.output_dir
    )

if __name__ == "__main__":
    main()
