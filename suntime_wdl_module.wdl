## WILDS WDL module for SJL (Solar Jetlag) suntime tile processing.
## Calculates sunrise/sunset time differences using one reference sunrise/ sunset per timezone for geographic tiles
## using NOAA solar calculator variables, as part of the suntime model pipeline.

version 1.0

task sjl_tiles {
  meta {
    author: "Caroline Nondin"
    email: "cnondin@fredhutch.org"
    description: "Calculate sunrise/sunset time differences for a geographic tile as part of the suntime model pipeline"
    url: "https://raw.githubusercontent.com/caronondin22/suntime_wdl/refs/heads/main/suntime_wdl_module.wdl"
    outputs: {
      matched_points: "RDS file containing points with sunrise/sunset difference values"
    }
    topic: "public_health_and_epidemiology"
    species: "human,eukaryote,prokaryote,virus"
    operation: "statistical_calculation"
    input_sample_required: "tile_path:accession:binary_format,border_points_path:accession:csv"
    input_sample_optional: "none"
    input_reference_required: "none"
    input_reference_optional: "none"
    output_sample: "matched_points:report:binary_format,missing_points:report:binary_format"
    output_reference: "none"
  }

  parameter_meta {
    tile_path: "Path to input tile .rds file"
    border_points_path: "Path to border points .csv file containing timezone boundary data"
    year: "Year for solar calculations (e.g. 2022)"
    matched_prefix: "Filename prefix for matched results output (default: 'matched_')"
    cpu_cores: "Number of CPU cores to use"
    memory_gb: "Memory allocation in GB"
  }

  input {
    File tile_path
    File border_points_path
    Int year
    String matched_prefix = "matched_"
    Int cpu_cores = 1
    Int memory_gb = 8
  }

  String tile_basename = basename(tile_path, ".rds")

  command <<<
    set -eo pipefail

    # Pull sjl_tiles script from GitHub
    # NOTE: For reproducibility in production workflows, replace the branch reference
    # (e.g., "refs/heads/main") with a specific commit hash (e.g., "abc1234...")
    wget -q "https://raw.githubusercontent.com/caronondin22/suntime_wdl/refs/heads/main/suntime_module.R" \
      -O sjl_tiles.R

    Rscript sjl_tiles.R \
      --tile_path "~{tile_path}" \
      --border_points_path "~{border_points_path}" \
      --year ~{year} \
      --matched_prefix "~{matched_prefix}" 
  >>>

  output {
    File matched_points = "~{matched_prefix}~{tile_basename}.rds"
  }

  runtime {
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
    docker: "getwilds/r-utils:0.1.0"
  }
}
