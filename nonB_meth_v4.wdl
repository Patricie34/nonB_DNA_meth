version 1.0
workflow nonb_methylation_complete_pipeline {
  input {
    File fasta                          
    Array[String] contig_filters = ["all"]  
    Array[String] motif_classes = ["G4"]    
    Int threads = 2
    String out_prefix = "hg002_nonB_meth"
    String extra_gfa_switches = ""
    Array[File] bedmethyl_files         
  }
  
  # Step 1: Split FASTA and get all contig names
  call split_fasta_task { input: fasta = fasta }
  
  # Step 2: Filter contigs based on contig_filters
  call filter_contigs {
    input: 
      all_contigs = split_fasta_task.contig_names,
      contig_filters = contig_filters
  }
  
  # Step 3: Process each contig for non-B DNA detection
  scatter (contig in filter_contigs.filtered_contigs) {
    call extract_contig_from_split {
      input: 
        fasta = fasta,
        contig = contig
    }
    call run_gfa {
      input:
        contig_fasta = extract_contig_from_split.contig_fasta,
        motif_list = if length(motif_classes) > 0 then motif_classes[0] else "G4",
        out_prefix = "${out_prefix}_${contig}",
        extra_gfa_switches = extra_gfa_switches,
        threads = threads
    }
  }
  
  # Step 4: Methylation analysis - process GFF files
  scatter (gff_idx in range(length(flatten(run_gfa.gff_files)))) {
    File gff_file = flatten(run_gfa.gff_files)[gff_idx]
    
    # Extract contig name from GFF filename - More robust extraction
    String contig_name_with_suffix = sub(sub(basename(gff_file), "\\.gff.*$", ""), "^~{out_prefix}_", "")
    String contig_name = sub(contig_name_with_suffix, "_GQ$", "")
    call gff_to_bed {
      input:
        gff = gff_file,
        contig = contig_name,
        out_prefix = out_prefix
    }
    
    # Process each bedMethyl file separately for this GFF
    scatter (bedmethyl_file in bedmethyl_files) {
      call filter_methylation_by_contig {
        input:
          bedmethyl_file = bedmethyl_file,
          contig = contig_name
      }
      
      call intersect_methyl_single {
        input:
          motif_bed = gff_to_bed.motif_bed,
          bedmethyl_file = filter_methylation_by_contig.filtered_bedmethyl
      }
    }
  }
  
  output {
    Array[String] processed_contigs = filter_contigs.filtered_contigs
    Array[File] gff_files = flatten(run_gfa.gff_files)
    Array[File] tsv_files = flatten(run_gfa.tsv_files)
    Array[File] bed_files = gff_to_bed.motif_bed
    Array[File] methylation_results = flatten(intersect_methyl_single.motif_methylation_summary)
    String workflow_summary = "Pipeline completed processing ${length(filter_contigs.filtered_contigs)} contigs"
  }
}

# -------------------------
# TASK: split_fasta_task (Simplified)
# Just extract contig names from FASTA
# -------------------------
task split_fasta_task {
  input {
    File fasta
  }
  command <<<
    set -eou pipefail
    mkdir -p split_output
    
    # Extract contig names directly from FASTA file
    grep "^>" "~{fasta}" | sed 's/^>//' | cut -d' ' -f1 > split_output/contig_names.txt
    
  >>>
  output {
    Array[String] contig_names = read_lines("split_output/contig_names.txt")
  }
  runtime {
    docker: "python:3.10-slim"
    cpu: 1
    memory: "2G"
  }
}

# -------------------------
# TASK: filter_contigs (SIMPLEST APPROACH)
# Filter contig names based on contig_filters
# -------------------------
task filter_contigs {
  input {
    Array[String] all_contigs
    Array[String] contig_filters
  }
  
  # Pre-process inputs using write_lines
  File contigs_file = write_lines(all_contigs)
  File filters_file = write_lines(contig_filters)
  
  command <<<
    set -eou pipefail
    
    # Convert line-separated files to comma-separated strings
    python3 - <<'PY1'
with open('~{contigs_file}', 'r') as f:
    contigs = [line.strip() for line in f if line.strip()]
contigs_str = ','.join(contigs)
with open('contigs_input.txt', 'w') as f:
    f.write(contigs_str)
PY1
    python3 - <<'PY2'  
with open('~{filters_file}', 'r') as f:
    filters = [line.strip() for line in f if line.strip()]
filters_str = ','.join(filters)
with open('filters_input.txt', 'w') as f:
    f.write(filters_str)
PY2
    # Main filtering logic
    python3 - <<'PY'
import sys
# Read inputs from files
with open('contigs_input.txt', 'r') as f:
    contigs_str = f.read().strip()
with open('filters_input.txt', 'r') as f:
    filters_str = f.read().strip()
# Parse contig names
if contigs_str and contigs_str != '':
    all_contigs = [x.strip() for x in contigs_str.split(',') if x.strip()]
else:
    all_contigs = []
# Parse filter patterns
if filters_str and filters_str != '':
    contig_filters = [x.strip() for x in filters_str.split(',') if x.strip()]
else:
    contig_filters = ['all']
# Apply filtering logic
if any(x.lower() == 'all' for x in contig_filters):
    result = all_contigs
else:
    result = []
    for filter_pattern in contig_filters:
        filter_pattern = filter_pattern.strip()
        if not filter_pattern:
            continue
        for contig in all_contigs:
            contig = contig.strip()
            if not contig:
                continue
            # Exact match
            if contig == filter_pattern:
                if contig not in result:
                    result.append(contig)
            # Prefix match
            elif contig.startswith(filter_pattern):
                if contig not in result:
                    result.append(contig)
# Output results
for contig in result:
    print(contig)
PY
  >>>
  output {
    Array[String] filtered_contigs = read_lines(stdout())
  }
  runtime {
    docker: "python:3.10-slim"
    cpu: 1
    memory: "2G"
  }
}

# -------------------------
# TASK: extract_contig_from_split (CORRECTED)
# Extract specific contig from split directory (simplified)
# -------------------------
task extract_contig_from_split {
  input {
    File fasta
    String contig
  }
  command <<<
    set -eou pipefail
    # Create a clear output directory structure
    OUTPUT_DIR="extracted_contigs"
    mkdir -p "${OUTPUT_DIR}"
    
    # Simple extraction using awk
    awk -v target_contig="~{contig}" '
    BEGIN { found = 0; target_header = ">" target_contig }
    /^>/ {
        if (found) exit
        if ($0 == target_header) found = 1
    }
    found { print }
    ' "~{fasta}" > "${OUTPUT_DIR}/~{contig}.fa"
    
    # Check if file was created and is not empty
    if [ ! -s "${OUTPUT_DIR}/~{contig}.fa" ]; then
      echo "ERROR: Contig ~{contig} not found or empty" >&2
      echo "Available contigs in FASTA:" >&2
      grep "^>" "~{fasta}" | head -10 >&2
      exit 1
    fi
    
    echo "Successfully created contig file" >&2
    ls -la "${OUTPUT_DIR}/" >&2
  >>>
  output {
    File contig_fasta = "extracted_contigs/~{contig}.fa"
  }
  runtime {
    docker: "python:3.10-slim"
    cpu: 1
    memory: "2G"
  }
}

# -------------------------
# TASK: run_gfa (FIXED)
# Runs gfa on a single-contig fasta and allows motif selection
# -------------------------
task run_gfa {
  input {
    File contig_fasta
    String motif_list        # comma-separated e.g. "G4" or "G4,MR" or "all"
    String out_prefix
    String extra_gfa_switches
    Int threads
  }
  command <<<
    set -eou pipefail
    mkdir -p out
    
    # DEBUG: Verify input file and show directory structure
    echo "=== INPUT VERIFICATION ===" >&2
    echo "Input contig_fasta: ~{contig_fasta}" >&2
    echo "File exists: $(test -f "~{contig_fasta}" && echo "YES" || echo "NO")" >&2
    if [ -f "~{contig_fasta}" ]; then
      echo "File size: $(wc -c "~{contig_fasta}" | cut -d' ' -f1)" >&2
      echo "First 2 lines:" >&2
      head -2 "~{contig_fasta}" >&2
    fi
    echo "Current directory contents:" >&2
    ls -la >&2
    echo "=== END INPUT VERIFICATION ===" >&2
    
    # Ensure the input file exists
    if [ ! -f "~{contig_fasta}" ]; then
      echo "FATAL ERROR: Input contig FASTA file not found: ~{contig_fasta}" >&2
      exit 1
    fi
    
    # Create necessary directories to prevent GFA path issues
    mkdir -p /mnt/disks/cromwell_root 2>/dev/null || true
    touch /mnt/disks/cromwell_root/stderr 2>/dev/null || true
    touch stderr 2>/dev/null || true
    
    # build skip flags: default skip everything, then remove skip for requested motifs
    ALL_FLAGS="-skipAPR -skipSTR -skipDR -skipMR -skipIR -skipGQ -skipZ -skipSlipped -skipCruciform -skipTriplex"
    ML="$(echo '~{motif_list}' | tr -d '[:space:]' | tr ',' ' ')"
    SKIP_FLAGS="${ALL_FLAGS}"
    
    if [ "${ML}" = "all" ] || [ "${ML}" = "ALL" ] || [ -z "${ML}" ]; then
      SKIP_FLAGS=""
    else
      for k in ${ML}; do
        kk=$(echo ${k} | tr '[:lower:]' '[:upper:]')
        case "${kk}" in
          "G4"|"GQ") SKIP_FLAGS=$(echo "${SKIP_FLAGS}" | sed 's/-skipGQ//g') ;;
          "IR") SKIP_FLAGS=$(echo "${SKIP_FLAGS}" | sed 's/-skipIR//g') ;;
          "MR") SKIP_FLAGS=$(echo "${SKIP_FLAGS}" | sed 's/-skipMR//g') ;;
          "DR") SKIP_FLAGS=$(echo "${SKIP_FLAGS}" | sed 's/-skipDR//g') ;;
          "APR") SKIP_FLAGS=$(echo "${SKIP_FLAGS}" | sed 's/-skipAPR//g') ;;
          "STR") SKIP_FLAGS=$(echo "${SKIP_FLAGS}" | sed 's/-skipSTR//g') ;;
          "Z") SKIP_FLAGS=$(echo "${SKIP_FLAGS}" | sed 's/-skipZ//g') ;;
          "SLIPPED") SKIP_FLAGS=$(echo "${SKIP_FLAGS}" | sed 's/-skipSlipped//g') ;;
          "CRUCIFORM") SKIP_FLAGS=$(echo "${SKIP_FLAGS}" | sed 's/-skipCruciform//g') ;;
          "TRIPLEX") SKIP_FLAGS=$(echo "${SKIP_FLAGS}" | sed 's/-skipTriplex//g') ;;
          *) echo "WARNING: unknown motif key ${k} (ignored)" >&2 ;;
        esac
      done
    fi
    
    # normalize SKIP_FLAGS
    SKIP_FLAGS=$(echo ${SKIP_FLAGS} | xargs)
    
    # call gfa with strict error handling
    OUTP="out/~{out_prefix}"
    echo "Running gfa: skip flags: ${SKIP_FLAGS}" >&2
    echo "Input file: ~{contig_fasta}" >&2
    echo "Output prefix: ${OUTP}" >&2
    
    # Run GFA with strict error handling - will fail the task if GFA fails
    /usr/local/bin/gfa -seq "~{contig_fasta}" -out "${OUTP}" ${SKIP_FLAGS} ~{extra_gfa_switches} 2>&1
    
    # Check what files were actually created
    echo "=== OUTPUT FILES CHECK ===" >&2
    if ls "out/~{out_prefix}"*.gff 1> /dev/null 2>&1; then
      echo "GFF files created:" >&2
      ls -la "out/~{out_prefix}"*.gff >&2
    else
      echo "No GFF files found - this may indicate GFA did not find motifs" >&2
    fi
    
    if ls "out/~{out_prefix}"*.tsv 1> /dev/null 2>&1; then
      echo "TSV files created:" >&2
      ls -la "out/~{out_prefix}"*.tsv >&2
    else
      echo "No TSV files found" >&2
    fi
    echo "=== END OUTPUT FILES CHECK ===" >&2
  >>>
  output {
    Array[File] gff_files = glob("out/~{out_prefix}*.gff")
    Array[File] tsv_files = glob("out/~{out_prefix}*.tsv")
  }
  runtime {
    docker: "patricie/nonb_gfa-image:v1.2" 
    cpu: 4
    memory: "16G"
  }
}

task gff_to_bed {
  input {
    File gff
    String contig
    String out_prefix
  }
  command <<<
    set -eou pipefail
    mkdir -p out
    
    # Debug: print what we received
    echo "DEBUG: gff = ~{gff}" >&2
    echo "DEBUG: contig = ~{contig}" >&2
    echo "DEBUG: out_prefix = ~{out_prefix}" >&2
    echo "DEBUG: gff basename = $(basename "~{gff}")" >&2
    
    # Check if file exists and show first few lines
    if [ -f "~{gff}" ]; then
      echo "DEBUG: First 5 lines of GFF file:" >&2
      head -5 "~{gff}" >&2
    else
      echo "ERROR: GFF file does not exist!" >&2
      exit 1
    fi
    
    # Write inputs to files to avoid shell escaping issues
    echo "~{gff}" > gff_input.txt
    echo "~{contig}" > contig_input.txt
    echo "~{out_prefix}" > out_prefix_input.txt
    
    python3 - <<'PY'
import sys
import re
import os
# Read inputs from files to avoid shell escaping issues
with open('gff_input.txt', 'r') as f:
    gff_file = f.read().strip()
with open('contig_input.txt', 'r') as f:
    contig = f.read().strip()
with open('out_prefix_input.txt', 'r') as f:
    out_prefix = f.read().strip()
print(f"DEBUG: Python received - gff_file: {gff_file}, contig: {contig}, out_prefix: {out_prefix}", file=sys.stderr)
# Validate inputs
if not os.path.exists(gff_file):
    print(f"ERROR: GFF file {gff_file} does not exist", file=sys.stderr)
    sys.exit(1)
gff_basename = os.path.basename(gff_file)
motif_class = gff_basename.split('_')[-1].replace('.gff', '')  # e.g., G4
out_file = f'out/{out_prefix}_{contig}_{motif_class}.bed'
print(f"DEBUG: Output file will be: {out_file}", file=sys.stderr)
line_count = 0
valid_lines = 0
try:
    with open(gff_file) as fh, open(out_file, 'w') as fo:
        for line in fh:
            line_count += 1
            line = line.rstrip('\n\r')
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
                
            parts = line.split('\t')
            
            # Must have exactly 9 fields
            if len(parts) != 9:
                print(f"WARNING: Line {line_count} has {len(parts)} fields, expected 9. Skipping: {line}", file=sys.stderr)
                continue
                
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            
            # Validate coordinates
            try:
                start_int = int(start)
                end_int = int(end)
                if start_int >= end_int:
                    print(f"WARNING: Invalid coordinates on line {line_count}: {start} >= {end}. Skipping: {line}", file=sys.stderr)
                    continue
            except ValueError:
                print(f"WARNING: Invalid coordinates on line {line_count}: {start}, {end}. Skipping: {line}", file=sys.stderr)
                continue
            
            # Extract motif name from ID or use feature type
            motif_name = ftype
            m = re.search(r'ID=([^;]+)', attrs)
            if m:
                motif_name = m.group(1)
            
            # Convert to 0-based coordinates for BED
            bed_start = start_int - 1
            bed_end = end_int
            
            # Use score if available, otherwise 0
            bed_score = score if score != '.' else '0'
            
            # Validate strand
            if strand not in ['+', '-', '.']:
                strand = '.'
            
            # Write BED format: chrom, start, end, name, score, strand
            fo.write(f"{seqid}\t{bed_start}\t{bed_end}\t{motif_name}\t{bed_score}\t{strand}\n")
            valid_lines += 1
    
    print(f"Processed {line_count} lines, wrote {valid_lines} valid entries to {out_file}", file=sys.stderr)
    print(out_file)
    
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)
    
PY
  >>>
  output {
    File motif_bed = read_string(stdout())
  }
  runtime {
    docker: "python:3.10-slim"
    cpu: 1
    memory: "2G"
  }
}

task filter_methylation_by_contig {
  input {
    File bedmethyl_file
    String contig
  }
  
  command <<<
    set -euo pipefail
    mkdir -p out
    
    # Extract basename of bedmethyl file without extension
    basename_bedmethyl=$(basename "~{bedmethyl_file}" .bed)
    
    echo "Filtering ~{bedmethyl_file} for contig ~{contig}" >&2
    
    # Correct: look for contig + TAB
    printf "^%s\t" "~{contig}" > pattern.txt
    
    # Filter and create descriptive filename
    grep -f pattern.txt "~{bedmethyl_file}" > out/${basename_bedmethyl}_~{contig}.bed || true
    
    echo "Filtered lines: $(wc -l < out/${basename_bedmethyl}_~{contig}.bed)" >&2
    echo "out/${basename_bedmethyl}_~{contig}.bed"
  >>>
  output {
    File filtered_bedmethyl = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:22.04"
    cpu: 1
    memory: "2G"
  }
}

task intersect_methyl_single {
  input {
    File motif_bed
    File bedmethyl_file
  }
  
  command <<<
    set -eou pipefail
    mkdir -p out
    
    # Extract basenames for descriptive naming
    basename_motif=$(basename "~{motif_bed}" .bed)
    basename_methyl=$(basename "~{bedmethyl_file}" .bed)
    
    # Check if files exist and have content
    if [ ! -s "~{motif_bed}" ] || [ ! -s "~{bedmethyl_file}" ]; then
        touch out/${basename_methyl}_intersect_${basename_motif}.bed
        echo "out/${basename_methyl}_intersect_${basename_motif}.bed"
        exit 0
    fi
    
    # Perform intersection - get bedMethyl entries that overlap with motifs
    # -a bedmethyl_file: the entries we want to output
    # -b motif_bed: the regions to intersect with  
    # -u: report each entry in A that overlaps with B (unique)
    bedtools intersect -a "~{bedmethyl_file}" -b "~{motif_bed}" -u > out/${basename_methyl}_intersect_${basename_motif}.bed 2>/dev/null || true
    
    echo "out/${basename_methyl}_intersect_${basename_motif}.bed"
  >>>
  output {
    File motif_methylation_summary = read_string(stdout())
  }
  runtime {
    docker: "patricie/bedtools_py-image:v1.2"
  }
}
