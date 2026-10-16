from pathlib import Path
import re
import pandas as pd
import argparse

def extract_pipeline_info(text):
    info = {}
    text = text.splitlines()
    
    # Dataset
    dataset_match = re.search(r"Working on\s*:\s*(.+)", text[0])
    info["Dataset"] = dataset_match.group(1).strip() if dataset_match else None

    # Loss Function
    loss_match = re.search(r"Estimation error function\s*:\s*(.+)", text[1])
    info["Loss Function"] = loss_match.group(1).strip() if loss_match else None
    
    # nbases
    nbases_match = re.search(r"nbases\s*:\s*(.+)", text[2])
    info["nbases"] = nbases_match.group(1).strip() if nbases_match else None
    
    # HPC node
    node_match = re.search(r"Running on master node:\s*(\S+)", text[4])
    info["HPC node"] = node_match.group(1).strip() if node_match else None
    
    # Total bases and reads used
    stats_match = re.search(r"(\d+)\s+total bases in\s+(\d+)\s+reads", text[9])
    if stats_match:
        info["Total bases used"] = int(stats_match.group(1))
        info["Total reads used"] = int(stats_match.group(2))
    else:
        info["Total bases used"] = None
        info["Total reads used"] = None
        
    return info

def parse_output_file(file_path):
    text = Path(file_path).read_text(encoding="utf-8")
    return extract_pipeline_info(text)

def parse_folder_to_csv(folder_path, output_csv="pipeline_summary.csv", file_pattern="*.out"):
    folder = Path(folder_path)
    results = []
    
    for file_path in folder.glob(file_pattern):
        if file_path.is_file():
            results.append(parse_output_file(file_path))
            
    # Convert to DataFrame and save to CSV
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"Successfully processed {len(results)} files and saved the summary to {output_csv}")
    return df

if __name__=="__main__":
    description = (
            "gather_logs.py : "
            "This script allow user to gather from dada2_workflow logs information about Dataset, Loss Function, nbases, HPC node, Total bases used, Total reads used."
        )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--datadir", "-d", required=True, help="Path to the directory with all output logs")
    parser.add_argument("--outfile", "-o", required=True, help="Path to the output csv file")
    args = parser.parse_args()
    parse_folder_to_csv(args.datadir, args.outfile)
