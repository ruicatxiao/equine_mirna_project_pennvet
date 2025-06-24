import os
import csv
import glob

def process_file(filepath):
    """
    Process a single file to extract the novel and known miRNA data.
    
    The file is assumed to be a tab-delimited text file with multiple sections.
    We look for two sections:
      - "novel miRNAs predicted by miRDeep2"
      - "mature miRBase miRNAs detected by miRDeep2"
    
    For the novel section, we extract:
      - ID from the "provisional id" column
      - Location from the "precursor coordinate" column
    
    For the known section, we extract:
      - ID from the "mature miRBase miRNA" column
      - Location from the "precursor coordinate" column
    
    The sample name is extracted from the file name, assuming it follows the pattern:
        result_SAMPLENAME.csv
    """
    # Extract sample name from file name "result_SAMPLENAME.csv"
    basename = os.path.basename(filepath)
    # Remove the "result_" prefix and ".csv" extension to get the sample name.
    sample = basename[len("result_"):-len(".csv")]
    
    # Define the section headers (converted to lower-case for easy matching)
    novel_section_header = "novel miRNAs predicted by miRDeep2".lower()
    known_section_header = "mature miRBase miRNAs detected by miRDeep2".lower()

    # This list will hold tuples: (sample, type, ID, location)
    results = []
    
    # Read the file content
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # Clean up lines: strip whitespace and ignore blank lines
    lines = [line.strip() for line in lines if line.strip() != ""]
    
    # Initialize variables for current section and header columns
    current_section = None  # will be "novel" or "known"
    header_indices = {}     # to store the column index for ID and location
    
    # Process each line in the file
    for i, line in enumerate(lines):
        # Use lower-case version for comparison
        lline = line.lower()
        
        # Check if the line marks the start of a section we care about.
        if novel_section_header in lline:
            current_section = "novel"
            header_indices = {}  # reset header mapping for this section
            continue
        elif known_section_header in lline:
            current_section = "known"
            header_indices = {}
            continue
        
        # If we have just switched into a section, then the next nonempty line should be the header row.
        if current_section is not None and not header_indices:
            headers = line.split("\t")
            if current_section == "novel":
                # Look for columns that contain "provisional id" and "precursor coordinate"
                for idx, h in enumerate(headers):
                    h_lower = h.lower()
                    if "provisional id" in h_lower:
                        header_indices["id"] = idx
                    if "precursor coordinate" in h_lower:
                        header_indices["loc"] = idx
            elif current_section == "known":
                # Look for columns that contain "mature mirbase mirna" and "precursor coordinate"
                for idx, h in enumerate(headers):
                    h_lower = h.lower()
                    if "mature mirbase mirna" in h_lower:
                        header_indices["id"] = idx
                    if "precursor coordinate" in h_lower:
                        header_indices["loc"] = idx
            continue  # after processing the header row, move to next line
        
        # If we are within a section and have identified the header, process the data rows.
        if current_section is not None and header_indices:
            # Split the row into fields using tab as the delimiter.
            data = line.split("\t")
            # Validate that there are enough columns in the row.
            if len(data) <= max(header_indices.values()):
                continue  # skip this line if it doesn't contain expected columns
            
            # Extract the ID and location values from the appropriate columns.
            id_val = data[header_indices["id"]].strip()
            loc_val = data[header_indices["loc"]].strip()
            # Skip empty or malformed rows.
            if not id_val or not loc_val:
                continue
            
            # Append the extracted row with the sample name and type.
            results.append((sample, current_section, id_val, loc_val))
    
    return results

def process_all_files(folder):
    """
    Process all files in the given folder that match the pattern 'result_*.csv'.
    For each file, a new report file will be generated named 'report_SAMPLENAME.txt'
    with the following columns:
      sample    type    ID    location
    """
    file_pattern = os.path.join(folder, "result_*.csv")
    all_files = glob.glob(file_pattern)
    all_reports = {}
    
    for file in all_files:
        results = process_file(file)
        if results:
            # Get the sample name from the first row (all rows in this file have the same sample)
            sample_name = results[0][0]
            # Define output file name (you can change the extension if desired)
            output_file = os.path.join(folder, f"report_{sample_name}.txt")
            with open(output_file, 'w', newline='', encoding='utf-8') as out_f:
                writer = csv.writer(out_f, delimiter="\t")
                # Write header line
                writer.writerow(["sample", "type", "ID", "location"])
                # Write each extracted row
                for row in results:
                    writer.writerow(row)
            all_reports[sample_name] = output_file
    return all_reports

if __name__ == "__main__":
    # Set the folder containing the mirdeep2 CSV files.
    folder = "."  # Change to your folder path if necessary
    reports = process_all_files(folder)
    print("Generated reports:")
    for sample, report in reports.items():
        print(f"{sample}: {report}")
