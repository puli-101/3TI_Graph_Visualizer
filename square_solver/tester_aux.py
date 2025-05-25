#!/usr/bin/env python3
import sys
import re
import csv

"""
    Auxiliary script to tester: Cleans up execution time and ressources output 
"""

def parse_time(time_str):
    """
    Convert a time string in h:mm:ss, m:ss, or even s format 
    (with optional fractional seconds) to seconds as a float.
    """
    parts = time_str.strip().split(':')
    try:
        if len(parts) == 3:
            # h:mm:ss(.ms)
            hours, minutes, seconds = parts
            return float(hours) * 3600 + float(minutes) * 60 + float(seconds)
        elif len(parts) == 2:
            # mm:ss(.ms)
            minutes, seconds = parts
            return float(minutes) * 60 + float(seconds)
        else:
            # If there's no colon, assume it's a seconds value.
            return float(time_str)
    except ValueError:
        # Fallback: If conversion fails, return None
        return None

def extract_metrics(block):
    """
    Extract the elapsed (wall clock) time and maximum resident set size
    from a block of /usr/bin/time -v output.
    Returns a tuple (elapsed_seconds, max_memory_kb) if both values are found.
    Otherwise, returns (None, None).
    """
    elapsed = None
    mem = None

    # Match the elapsed time line.
    # Example line:
    # "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.71"
    elapsed_match = re.search(r"Elapsed \(wall clock\) time.*:\s*([0-9:.]+)", block)
    if elapsed_match:
        time_str = elapsed_match.group(1)
        elapsed = parse_time(time_str)

    # Match the maximum resident set size.
    # Example line:
    # "Maximum resident set size (kbytes): 197832"
    mem_match = re.search(r"Maximum resident set size.*:\s*(\d+)", block)
    if mem_match:
        mem = int(mem_match.group(1))

    return elapsed, mem

def main():
    # Check command-line arguments
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: python3 parse_time_output.py <input_file>\n")
        sys.exit(1)
        
    input_file = sys.argv[1]
    try:
        with open(input_file, 'r') as f:
            content = f.read()
    except Exception as e:
        sys.stderr.write(f"Error opening file {input_file}: {e}\n")
        sys.exit(1)

    # If the file contains multiple /usr/bin/time -v outputs, they often start with
    # "Command being timed:" so we can split the content into blocks.
    blocks = re.split(r"\n(?=Command being timed:)", content)
    results = []
    for block in blocks:
        elapsed, mem = extract_metrics(block)
        if elapsed is not None and mem is not None:
            results.append((elapsed, mem))

    # Write CSV table to stdout
    writer = csv.writer(sys.stdout)
    #writer.writerow(["Execution Time (seconds)", "Max Memory Usage (kB)"])
    for elapsed, mem in results:
        writer.writerow([elapsed, mem])

if __name__ == "__main__":
    main()
