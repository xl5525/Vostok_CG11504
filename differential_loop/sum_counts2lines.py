import pandas as pd
import re
import sys

def process_file(input_file):
    # Extract the base filename without extension
    out_file = re.search(r'(.+?).txt.csv', input_file).group(1)
    
    merged_rows = []

    # Read the file line by line
    with open(input_file, 'r') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 2):
            # Ensure that we have an even number of lines
            if i + 1 < len(lines):
                # Merge two consecutive lines
                merged_row = lines[i].strip().split() + lines[i + 1].strip().split()
                merged_rows.append(merged_row)

    # Convert merged rows to DataFrame
    merged_df = pd.DataFrame(merged_rows)

    # Ensure the 14th column (index 13) is numeric, filling non-numeric values with 0
    merged_df[13] = pd.to_numeric(merged_df[13], errors='coerce').fillna(0)

    # Group by the 7th column (loop_id) and sum the 14th column
    df_grouped = merged_df.groupby(6)[13].sum().reset_index()

    # Rename the columns
    df_grouped.columns = ['loop_id', out_file]

    # Write the output to a CSV file
    df_grouped.to_csv(f'{out_file}_summed.csv', sep='\t', index=False)

if __name__ == "__main__":
    # Read the input file from command line argument
    input_file = sys.argv[1]

    # Process the file
    process_file(input_file)
