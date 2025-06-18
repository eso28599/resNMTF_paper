 # code obtained from chatgpt
import pandas as pd
import os

def convert_excel_to_csvs(filepath):
    # Load the Excel file
    xls = pd.ExcelFile(filepath)

    # Get the directory and filename base
    directory = os.path.dirname(filepath)
    base_name = os.path.splitext(os.path.basename(filepath))[0]

    # Iterate through each sheet and save as a CSV
    for sheet_name in xls.sheet_names:
        df = xls.parse(sheet_name)
        # Create a safe filename using the base name and sheet name
        safe_sheet_name = sheet_name.replace(" ", "_").replace("/", "_")
        csv_filename = f"{base_name}_{safe_sheet_name}.csv"
        csv_path = os.path.join(directory, csv_filename)
        df.to_csv(csv_path, index=False)
        print(f"Saved sheet '{sheet_name}' to '{csv_path}'")


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python excel_to_csvs.py <excel_file_path>")
    else:
        convert_excel_to_csvs(sys.argv[1])