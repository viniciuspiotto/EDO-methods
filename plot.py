import matplotlib.pyplot as plt
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
            method = "Exato"

        block = []
        i += 1
        while i < len(lines):
            l = lines[i].strip()
            if not l or (not (l[0].isdigit() or l[0] in ['-', '.'])):
                break
            block.append(l)
            i += 1

        data[method] = block

    for method, block in data.items():
        xs, ys = zip(*[map(float, line.split()) for line in block[:-1]])
        error = float(block[-1])
        data[method] = {'x': list(xs), 'y': list(ys), 'erro': error}

    return data

def plotIndividualGraph(x, y, exactX, exactY, label, exactLabel='Exato'):
    plt.figure(figsize=(10,6))
    plt.plot(exactX, exactY, label=exactLabel, marker='o')
    plt.plot(x, y, label=label, marker='o', linestyle='--')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Tensão do Capacitor (V)')
    plt.title(label)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    os.makedirs("graficos", exist_ok=True)
    plt.savefig(f'graficos/{label.lower().replace(" ", "")}.png')
    plt.close()

def plotGraphs(data):
    os.makedirs("graficos", exist_ok=True)
    exact = data.pop("Exato")

    plt.figure(figsize=(10,6))
    names = list(data.keys())
    errors = [v['erro'] for v in data.values()]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
    bars = plt.bar(names, errors, color=colors[:len(names)], width=0.6)
    plt.yscale('log')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Erro RMS (escala log)')
    plt.title('Erro RMS dos Métodos')
    plt.subplots_adjust(bottom=0.25)
    ymin = plt.ylim()[0]
    for bar, error in zip(bars, errors):
        x = bar.get_x() + bar.get_width() / 2
        plt.text(x + 0.1, ymin * 0.8, f"{abs(error):.2f}", ha='center', va='top', fontsize=9, rotation=45)
    plt.tight_layout()
    plt.savefig('graficos/erros.png')
    plt.close()

    stableMethods = ["Euler Implicito", "BDF2", "Adams Moulton 2", "Trapezio Implicito"]
    plt.figure(figsize=(10,6))
    plt.plot(exact['x'], exact['y'], label='Exato', marker='o')
    for method in stableMethods:
        if method in data:
            plt.plot(data[method]['x'], data[method]['y'], label=method, marker='o', linestyle='--')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Tensão do Capacitor (V)')
    plt.title('Métodos Estáveis')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('graficos/metodosEstaveis.png')
    plt.close()

    for method in ["Euler Explicito", "Adams Bashford 2"]:
        if method in data:
            plotIndividualGraph(data[method]['x'], data[method]['y'], exact['x'], exact['y'], method)

data = readFile("EDO_Descarga_de_Um_Capacitor_de_Circuitos_de_RC.txt")
plotGraphs(data)
