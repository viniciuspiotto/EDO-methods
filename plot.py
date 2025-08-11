import matplotlib.pyplot as plt

def ler_arquivo_e_plotar(nome_arquivo):
    dados = {}
    metodo_atual = None

    with open(nome_arquivo, 'r') as f:
        for linha in f:
            linha = linha.strip()
            if not linha:
                continue
            if not linha[0].isdigit() and not linha[0] == '-' and not linha[0] == '.':
                metodo_atual = linha
                dados[metodo_atual] = {'x': [], 'y': []}
            else:
                x_str, y_str = linha.split()
                x = float(x_str)
                y = float(y_str)
                dados[metodo_atual]['x'].append(x)
                dados[metodo_atual]['y'].append(y)

    plt.figure(figsize=(12, 8))

    for metodo, valores in dados.items():
        plt.plot(valores['x'], valores['y'], label=metodo)

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Soluções numéricas dos métodos')
    plt.legend()
    plt.grid(True)
    plt.show()

ler_arquivo_e_plotar('EDO_Descarga_de_Um_Capacitor_de_Circuitos_de_RC.txt')
