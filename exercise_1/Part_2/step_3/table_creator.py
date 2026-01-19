import re
import os
from pathlib import Path

def parse_filename(filename):
    """Extract topology and grid size from filename."""
    pattern = r'output_t(\d+)x(\d+)_g(\d+)x(\d+)_i(\d+)'
    match = re.match(pattern, filename)
    if match:
        return {
            'topo_x': int(match.group(1)),
            'topo_y': int(match.group(2)),
            'grid_x': int(match.group(3)),
            'grid_y': int(match.group(4)),
            'max_iter': int(match.group(5))
        }
    return None

def parse_logfile(filepath):
    """Extract wtime and iterations from log file."""
    wtimes = []
    iterations = []
    
    with open(filepath, 'r') as f:
        for line in f:
            # Look for iteration count
            if 'Number of iterations' in line:
                match = re.search(r'Number of iterations\s*:\s*(\d+)', line)
                if match:
                    iterations.append(int(match.group(1)))
            
            # Look for wtime
            if 'Elapsed Wtime' in line:
                match = re.search(r'Elapsed Wtime\s+([\d.]+)', line)
                if match:
                    wtimes.append(float(match.group(1)))
    
    # Return max wtime and check iterations are consistent
    wtime = max(wtimes) if wtimes else None
    iter_count = iterations[0] if iterations and all(i == iterations[0] for i in iterations) else None
    
    return wtime, iter_count

def create_latex_table(data, output_file='results_table.tex'):
    """Create a LaTeX table from the parsed data."""
    
    # Sort data by topology, then grid size, then max iterations
    sorted_data = sorted(data, key=lambda x: (x['topo_x'], x['topo_y'], 
                                               x['grid_x'], x['grid_y'],
                                               x['max_iter']))
    
    latex = []
    latex.append(r'\begin{table}[htbp]')
    latex.append(r'\centering')
    latex.append(r'\begin{tabular}{cccccc}')
    latex.append(r'\hline')
    latex.append(r'Topology & Grid Size & Max Iter & Actual Iter & Wall Time (s) & Time/Iter (ms) \\')
    latex.append(r'\hline')
    
    for row in sorted_data:
        topo = f"{row['topo_x']}$\\times${row['topo_y']}"
        grid = f"{row['grid_x']}$\\times${row['grid_y']}"
        max_iter = row['max_iter']
        actual_iter = row['iterations'] if row['iterations'] else 'N/A'
        wtime = f"{row['wtime']:.4f}" if row['wtime'] else 'N/A'
        
        # Calculate time per iteration if both are available
        if row['wtime'] and row['iterations']:
            time_per_iter = (row['wtime'] / row['iterations']) * 1000  # convert to ms
            time_per_iter_str = f"{time_per_iter:.2f}"
        else:
            time_per_iter_str = 'N/A'
        
        latex.append(f"{topo} & {grid} & {max_iter} & {actual_iter} & {wtime} & {time_per_iter_str} \\\\")
    
    latex.append(r'\hline')
    latex.append(r'\end{tabular}')
    latex.append(r'\caption{Poisson solver performance results}')
    latex.append(r'\label{tab:results}')
    latex.append(r'\end{table}')
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(latex))
    
    print(f"\nLaTeX table written to {output_file}")
    return '\n'.join(latex)

def main():
    # Directory containing the log files
    log_dir = './'  # Adjust this path as needed
    
    data = []
    
    # Process all log files
    for filename in sorted(os.listdir(log_dir)):
        if filename.startswith('output_t') and filename.endswith('.log'):
            # Parse filename
            file_params = parse_filename(filename)
            if not file_params:
                print(f"Warning: Could not parse filename {filename}")
                continue
            
            # Parse log file content
            filepath = os.path.join(log_dir, filename)
            wtime, iterations = parse_logfile(filepath)
            
            # Combine data
            file_params['wtime'] = wtime
            file_params['iterations'] = iterations
            file_params['filename'] = filename
            
            data.append(file_params)
            print(f"Processed: {filename:<50} wtime: {wtime:>8.4f}s, iterations: {iterations}")
    
    if not data:
        print("No log files found!")
        return
    
    # Create LaTeX table
    latex_table = create_latex_table(data)
    print("\n" + "="*80)
    print("Generated LaTeX table:")
    print("="*80)
    print(latex_table)

if __name__ == "__main__":
    main()