import re
import os
import numpy as np

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
            if 'Number of iterations' in line:
                match = re.search(r'Number of iterations\s*:\s*(\d+)', line)
                if match:
                    iterations.append(int(match.group(1)))
            
            if 'Elapsed Wtime' in line:
                match = re.search(r'Elapsed Wtime\s+([\d.]+)', line)
                if match:
                    wtimes.append(float(match.group(1)))
    
    wtime = max(wtimes) if wtimes else None
    iter_count = iterations[0] if iterations and all(i == iterations[0] for i in iterations) else None
    
    return wtime, iter_count

def main():
    log_dir = './'
    
    data = []
    
    # Process all log files
    for filename in sorted(os.listdir(log_dir)):
        if filename.startswith('output_t') and filename.endswith('.log'):
            file_params = parse_filename(filename)
            if not file_params:
                continue
            
            filepath = os.path.join(log_dir, filename)
            wtime, iterations = parse_logfile(filepath)
            
            if wtime is not None and iterations is not None:
                data.append({
                    'topology': f"{file_params['topo_x']}x{file_params['topo_y']}",
                    'grid_size': f"{file_params['grid_x']}x{file_params['grid_y']}",
                    'grid_x': file_params['grid_x'],
                    'grid_y': file_params['grid_y'],
                    'iterations': iterations,
                    'wall_time': wtime
                })
    
    # Sort by topology, grid size, then iterations
    data.sort(key=lambda x: (x['topology'], x['grid_x'], x['iterations']))
    
    # Print header
    print(f"{'Topology':<10} {'Grid Size':<12} {'Iterations':<12} {'Wall Time (s)':<15}")
    print("-" * 50)
    
    # Print data
    for d in data:
        print(f"{d['topology']:<10} {d['grid_size']:<12} {d['iterations']:<12} {d['wall_time']:<15.6f}")
    
    # Save to CSV
    with open('performance_data.csv', 'w') as f:
        f.write("topology,grid_size,iterations,wall_time\n")
        for d in data:
            f.write(f"{d['topology']},{d['grid_size']},{d['iterations']},{d['wall_time']:.6f}\n")
    
    print("\nData saved to 'performance_data.csv'")

if __name__ == "__main__":
    main()