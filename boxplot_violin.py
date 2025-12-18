#!/usr/bin/env python3
"""
Methylation Analysis at Non-B DNA Motifs

This script analyzes methylation levels at various non-B DNA motifs by processing
intersected BED files and generating comparative visualizations.

Author: Your Name
License: MIT
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import argparse
import json
import sys

def create_methylation_analysis(input_files, output_prefix="methylation_analysis", 
                              chromosome="chrX", haplotype="Maternal", output_dir="./plots"):
    """
    Analyze methylation levels at non-B DNA motifs
    
    Parameters:
    input_files (dict or list): Dictionary with descriptive keys and file paths as values,
                               or list of file paths
    output_prefix (str): Prefix for output files
    chromosome (str): Chromosome name
    haplotype (str): Haplotype name (Maternal/Paternal)
    output_dir (str): Directory to save output files
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Define column names for BEDMethyl format
    columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 
               'thickEnd', 'itemRgb', 'coverage', 'methylation_pct', 'methylated', 
               'unmethylated', 'context1', 'context2', 'context3', 'context4', 'context5']
    
    # Motif mapping - maps file identifiers to display names
    motif_mapping = {
        'G4': 'G-Quadruplex (G4)',
        'GQ': 'G-Quadruplex (GQ)',
        'IR': 'Inverted Repeat (IR)',
        'MR': 'Mirror Repeat (MR)',
        'DR': 'Direct Repeat (DR)',
        'APR': 'APR Repeat (APR)',
        'STR': 'Short Tandem Repeat (STR)',
        'Z': 'Z-DNA (Z)',
        'SLIPPED': 'Slipped Structure (Slipped)',
        'CRUCIFORM': 'Cruciform (Cruciform)',
        'TRIPLEX': 'Triplex (Triplex)'
    }
    
    # Dictionary to store dataframes
    dataframes = {}
    
    # Handle both dict and list input formats
    if isinstance(input_files, dict):
        files_to_process = input_files.items()
    elif isinstance(input_files, list):
        # Create keys from filenames
        files_to_process = [(os.path.basename(f).replace('.bed', ''), f) for f in input_files]
    else:
        raise ValueError("input_files must be a dictionary or list")
    
    # Load the data
    for key, file_path in files_to_process:
        if file_path and os.path.exists(file_path):
            print(f"Loading {key} data...")
            try:
                df = pd.read_csv(file_path, sep='\t', header=None, names=columns)
                
                # Extract motif and modification types from key or filename
                motif_type, mod_type = extract_motif_mod_types(key, file_path, motif_mapping)
                
                df['motif_type'] = motif_type
                df['modification_type'] = mod_type
                
                dataframes[key] = df
            except Exception as e:
                print(f"Error loading {key}: {e}")
                dataframes[key] = pd.DataFrame()
        else:
            print(f"File not found or not specified for {key}")
            dataframes[key] = pd.DataFrame()
    
    # Check if we have any data
    if not dataframes or all(len(df) == 0 for df in dataframes.values()):
        print("No files found or no data loaded successfully.")
        return
    
    # Combine available dataframes
    df_list = [df for df in dataframes.values() if len(df) > 0]
    if not df_list:
        print("No data loaded successfully.")
        return
    
    df_combined = pd.concat(df_list, ignore_index=True)
    
    # Create combined categories for plotting
    df_combined['motif_mod_type'] = df_combined['motif_type'] + '\n' + df_combined['modification_type']
    
    # Get unique combinations and create color palette
    unique_combinations = df_combined['motif_mod_type'].unique()
    
    # Generate colors dynamically based on number of combinations
    colors = generate_color_palette(len(unique_combinations))
    color_map = dict(zip(unique_combinations, colors))
    
    # Create comparison boxplot
    plt.figure(figsize=(max(10, len(unique_combinations) * 2), 8))
    sns.boxplot(data=df_combined, x='motif_mod_type', y='methylation_pct', 
                palette=[color_map[combo] for combo in unique_combinations], width=0.6)
    sns.stripplot(data=df_combined, x='motif_mod_type', y='methylation_pct', 
                  color='black', alpha=0.4, size=3, jitter=0.2)
    plt.ylabel('Methylation Percentage (%)')
    plt.title(f'Methylation Levels at Non-B DNA Motifs\n{haplotype} Haplotype - {chromosome}')
    plt.grid(True, alpha=0.3)
    plt.xticks(rotation=45, ha='right')
    
    # Add sample sizes to the plot
    sample_sizes = df_combined['motif_mod_type'].value_counts()
    for i, label in enumerate(unique_combinations):
        n = sample_sizes.get(label, 0)
        plt.text(i, plt.ylim()[1] * 0.95, f'n={n}', ha='center', va='top', 
                 fontsize=8, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    output_filename = os.path.join(output_dir, f'{output_prefix}_comparison_{chromosome}.png')
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save summary statistics to CSV
    summary_data = []
    for key, df in dataframes.items():
        if len(df) > 0:
            motif_type = df['motif_type'].iloc[0]
            mod_type = df['modification_type'].iloc[0]
            summary_data.append({
                'File_Key': key,
                'Motif_Type': motif_type,
                'Modification_Type': mod_type,
                'Count': len(df),
                'Mean_Methylation': df['methylation_pct'].mean(),
                'Median_Methylation': df['methylation_pct'].median(),
                'Std_Deviation': df['methylation_pct'].std(),
                'Min_Methylation': df['methylation_pct'].min(),
                'Max_Methylation': df['methylation_pct'].max()
            })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        csv_filename = os.path.join(output_dir, f'{output_prefix}_summary_{chromosome}.csv')
        summary_df.to_csv(csv_filename, index=False)
    
    # Create additional visualization: violin plot
    plt.figure(figsize=(max(10, len(unique_combinations) * 2), 8))
    sns.violinplot(data=df_combined, x='motif_mod_type', y='methylation_pct', 
                   palette=[color_map[combo] for combo in unique_combinations])
    plt.ylabel('Methylation Percentage (%)')
    plt.title(f'Methylation Distribution at Non-B DNA Motifs\n{haplotype} Haplotype - {chromosome}')
    plt.grid(True, alpha=0.3)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    violin_filename = os.path.join(output_dir, f'{output_prefix}_violin_{chromosome}.png')
    plt.savefig(violin_filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create separate plots for different modification types if we have multiple motifs
    available_motifs = set()
    available_mods = set()
    
    for key, df in dataframes.items():
        if len(df) > 0:
            motif_type = df['motif_type'].iloc[0]
            mod_type = df['modification_type'].iloc[0]
            available_motifs.add(motif_type)
            available_mods.add(mod_type)
    
    if len(available_motifs) > 1:
        mod_types_available = list(available_mods)
        if len(mod_types_available) > 0:
            fig, axes = plt.subplots(1, len(mod_types_available), figsize=(6*len(mod_types_available), 6))
            if len(mod_types_available) == 1:
                axes = [axes]
            
            for idx, mod_type in enumerate(mod_types_available):
                # Get dataframes for this modification type
                mod_dfs = {key: df for key, df in dataframes.items() 
                          if len(df) > 0 and df['modification_type'].iloc[0] == mod_type}
                
                if mod_dfs:
                    df_mod = pd.concat(list(mod_dfs.values()), ignore_index=True)
                    motif_types = df_mod['motif_type'].unique()
                    
                    # Generate colors for motifs
                    motif_colors = generate_color_palette(len(motif_types))
                    
                    sns.boxplot(data=df_mod, x='motif_type', y='methylation_pct', 
                                palette=motif_colors, width=0.5, ax=axes[idx])
                    sns.stripplot(data=df_mod, x='motif_type', y='methylation_pct', 
                                  color='black', alpha=0.4, size=4, jitter=0.2, ax=axes[idx])
                    axes[idx].set_ylabel('Methylation Percentage (%)')
                    axes[idx].set_title(f'{mod_type} Methylation Levels\n{haplotype} Haplotype - {chromosome}')
                    axes[idx].grid(True, alpha=0.3)
                    axes[idx].tick_params(axis='x', rotation=45)
                    
                    # Add sample sizes
                    sample_sizes_mod = df_mod['motif_type'].value_counts()
                    for i, motif in enumerate(motif_types):
                        n = sample_sizes_mod.get(motif, 0)
                        axes[idx].text(i, axes[idx].get_ylim()[1] * 0.95, f'n={n}', ha='center', va='top', 
                                       fontsize=10, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                else:
                    axes[idx].set_title(f'{mod_type} Data Not Available')
                    axes[idx].text(0.5, 0.5, 'No data', ha='center', va='center', transform=axes[idx].transAxes)
            
            plt.tight_layout()
            separate_filename = os.path.join(output_dir, f'{output_prefix}_comparison_separate_{chromosome}.png')
            plt.savefig(separate_filename, dpi=300, bbox_inches='tight')
            plt.close()
    
    # Perform statistical comparisons and save to text file
    stats_filename = os.path.join(output_dir, f'{output_prefix}_statistics_{chromosome}.txt')
    with open(stats_filename, 'w') as f:
        f.write("=== STATISTICAL COMPARISONS ===\n\n")
        
        # Compare between motifs for each modification type
        for mod_type in available_mods:
            motif_dfs = {}
            for key, df in dataframes.items():
                if len(df) > 0 and df['modification_type'].iloc[0] == mod_type:
                    motif_type = df['motif_type'].iloc[0]
                    motif_dfs[motif_type] = df
            
            motif_types = list(motif_dfs.keys())
            if len(motif_types) >= 2:
                # Compare first two motifs
                df_list = list(motif_dfs.values())
                if len(df_list) >= 2:
                    df1, df2 = df_list[0], df_list[1]
                    motif1, motif2 = list(motif_dfs.keys())[0], list(motif_dfs.keys())[1]
                    try:
                        t_stat, p_value = stats.ttest_ind(df1['methylation_pct'], df2['methylation_pct'])
                        f.write(f"{motif1} vs {motif2} {mod_type} methylation: p = {p_value:.2e}\n")
                    except Exception as e:
                        f.write(f"Could not compare {motif1} vs {motif2} {mod_type}: {e}\n")
        
        # Compare modifications within each motif type
        for motif_type in available_motifs:
            mod_dfs = {}
            for key, df in dataframes.items():
                if len(df) > 0 and df['motif_type'].iloc[0] == motif_type:
                    mod_type = df['modification_type'].iloc[0]
                    mod_dfs[mod_type] = df
            
            mod_types = list(mod_dfs.keys())
            if len(mod_types) >= 2:
                # Compare first two modifications
                df_list = list(mod_dfs.values())
                if len(df_list) >= 2:
                    df1, df2 = df_list[0], df_list[1]
                    mod1, mod2 = list(mod_dfs.keys())[0], list(mod_dfs.keys())[1]
                    try:
                        t_stat, p_value = stats.ttest_ind(df1['methylation_pct'], df2['methylation_pct'])
                        f.write(f"{motif_type} {mod1} vs {mod2}: p = {p_value:.2e}\n")
                    except Exception as e:
                        f.write(f"Could not compare {motif_type} {mod1} vs {mod2}: {e}\n")
    
    print(f"Analysis complete! All output files saved to: {output_dir}")
    print(f"Available motifs: {', '.join(sorted(available_motifs))}")
    print(f"Available modifications: {', '.join(sorted(available_mods))}")

def extract_motif_mod_types(key, file_path, motif_mapping):
    """
    Extract motif and modification types from key or filename
    """
    # Convert to uppercase for matching
    key_upper = key.upper()
    filename_upper = os.path.basename(file_path).upper()
    combined_text = key_upper + " " + filename_upper
    
    # Extract modification type (5mC or 5hmC)
    mod_type = "5mC"
    if "5HMC" in combined_text or "HMC" in combined_text:
        mod_type = "5hmC"
    elif "5MC" in combined_text:
        mod_type = "5mC"
    
    # Extract motif type
    motif_type = "Unknown Motif"
    
    # Check for motif patterns in the combined text
    for motif_key, display_name in motif_mapping.items():
        if motif_key.upper() in combined_text:
            motif_type = display_name
            break
    
    # If still unknown, try to extract from key/filename
    if motif_type == "Unknown Motif":
        # Look for common motif patterns
        if "GQ" in combined_text:
            motif_type = "G-Quadruplex (GQ)"
        elif "Z" in combined_text and "Z" not in "".join(list(motif_mapping.keys())):
            motif_type = "Z-DNA (Z)"
        elif "STR" in combined_text:
            motif_type = "Short Tandem Repeat (STR)"
        else:
            # Use the key as motif name
            motif_type = key.replace("_5mC", "").replace("_5hmC", "").replace("_", " ").title()
    
    return motif_type, mod_type

def generate_color_palette(n_colors):
    """
    Generate a color palette with n distinct colors
    """
    if n_colors <= 0:
        return []
    
    # Predefined color schemes for common cases
    predefined_palettes = {
        1: ['#66c2a5'],
        2: ['#66c2a5', '#fc8d62'],
        3: ['#66c2a5', '#fc8d62', '#8da0cb'],
        4: ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3'],
        5: ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854'],
        6: ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f'],
        7: ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494'],
        8: ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3']
    }
    
    if n_colors in predefined_palettes:
        return predefined_palettes[n_colors]
    
    # For larger numbers, generate colors using matplotlib
    import matplotlib.cm as cm
    import numpy as np
    
    colors = cm.tab10(np.linspace(0, 1, min(n_colors, 10))).tolist()
    if n_colors > 10:
        # Add more colors if needed
        additional_colors = cm.Set3(np.linspace(0, 1, n_colors - 10)).tolist()
        colors.extend(additional_colors)
    
    # Convert to hex colors
    hex_colors = [f'#{int(c[0]*255):02x}{int(c[1]*255):02x}{int(c[2]*255):02x}' for c in colors[:n_colors]]
    
    return hex_colors

def main():
    parser = argparse.ArgumentParser(description='Analyze methylation levels at non-B DNA motifs')
    parser.add_argument('--input-files', nargs='+', required=True,
                       help='Input BED files with methylation data')
    parser.add_argument('--output-prefix', default='methylation_analysis',
                       help='Prefix for output files')
    parser.add_argument('--chromosome', default='chrX',
                       help='Chromosome name')
    parser.add_argument('--haplotype', default='Maternal',
                       help='Haplotype name (Maternal/Paternal)')
    parser.add_argument('--output-dir', default='./plots',
                       help='Directory to save output files')
    
    args = parser.parse_args()
    
    # Create input files dictionary from command line arguments
    input_files = args.input_files
    
    create_methylation_analysis(
        input_files=input_files,
        output_prefix=args.output_prefix,
        chromosome=args.chromosome,
        haplotype=args.haplotype,
        output_dir=args.output_dir
    )

if __name__ == "__main__":
    main()
