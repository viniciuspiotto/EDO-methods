#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import os

def readFile(fileName):
    data = {}
    with open(fileName, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        
        if not (line[0].isdigit() or line[0] in ['-', '.']):
            method = line
        else:
            method = "Exact"

        while i < len(lines) and lines[i].strip() != "Results:":
            i += 1
        i += 1

        block = []
        while i < len(lines):
            l = lines[i].strip()
            if not l or not (l[0].isdigit() or l[0] == '-'):
                break
            block.append(l)
            i += 1

        xs, ys = zip(*[map(float, line.split()) for line in block])
        data[method] = {'x': list(xs), 'y': list(ys)}

    if "Exact" not in data:
        x_tmp = [
            0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25,
            2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0
        ]
        y_tmp = [
            0, 4.69951014480389, -8.66145140695345, -5.41048537162463, 8.60309100643824,
            5.40569485822865, -8.60348423572376, -5.40572713645401, 8.60348158616568,
            5.40572691896504, -8.60348160401826, -5.40572692043047, 8.60348160389797,
            5.40572692042060, -8.60348160389878, -5.40572692042066, 8.60348160389877,
            5.40572692042066, -8.60348160389877, -5.40572692042066, 8.60348160389877
        ]
        data["Exact"] = {'x': x_tmp, 'y': y_tmp, 'error': 0.0}

    return data

def plotGraph(data):
    os.makedirs("plots", exist_ok=True)
    exact = data.pop("Exact")

    stableMethods = ["Implicit Euler", "BDF2", "Adams-Moulton 2", "Implicit Trapezoidal"]
    plt.figure(figsize=(10,6))
    plt.plot(exact['x'], exact['y'], label='Exact', marker='o')
    for method in stableMethods:
        if method in data:
            plt.plot(data[method]['x'], data[method]['y'], label=method, marker='o', linestyle='--')
    plt.xlabel('Time (s)')
    plt.ylabel('Capacitor Voltage (V)')
    plt.title('Stable Methods')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('plots/stablemethods.png')
    plt.close()

if (len(sys.argv) < 2):
    print("Usage: ./plot.py <output_file>")
    exit(1)
data = readFile(sys.argv[1])
plotGraph(data)