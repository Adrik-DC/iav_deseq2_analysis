# ------------------------------------------
# Definir los parametros
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: Guardar los parametros de entrada
# Entrada: input, output, columna_gene, columna_lfc, columna_padj
# Salida: args
# ------------------------------------------

# import librerias
import os
import argparse


def parsear_argumentos():
    """Parsea los argumentos de la línea de comandos.

    Returns:
        argparse.Namespace: Argumentos parseados con atributos
            `input`, `output`
    """

    parser = argparse.ArgumentParser(description=" Analisis de IAV")

    parser.add_argument("-i", "--input", help="Ruta al archivo fasta de entrada")

    # parser.add_argument("-o", "--output", help="Ruta al archivo fasta de salida")

    parser.add_argument(
        "--columna-gene",
        type=int,
        default=1,
        help="Índice de la columna que contiene el nombre del gen",
    )

    parser.add_argument(
        "--columna-lfc",
        type=int,
        default=3,
        help="Índice de la columna que contiene el log2FoldChange",
    )

    parser.add_argument(
        "--columna-padj",
        type=int,
        default=7,
        help="Índice de la columna que contiene el padj",
    )

    args = parser.parse_args()

    return args


# ------------------------------------------
# Extraer los parametros de las columnas del archivo de entrada
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: Guardar la informacion necesaria de las columnas del archivo de entrada
# Entrada: input, columna_gene, columna_lfc, columna_padj
# Salida: interactions
# ------------------------------------------


argumentos = parsear_argumentos()
filename = argumentos.input
gene_col = argumentos.columna_gene
lfc_col = argumentos.columna_lfc
padj_col = argumentos.columna_padj
print(argumentos)


def load_deseq2_results(filename, gene_col, lfc_col, padj_col):
    interactions = []
    if not os.path.exists(filename):
        print("Error: archivo no encontrado")
        exit(1)
    index = max(gene_col, lfc_col, padj_col)
    with open(filename) as f:
        f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0].islower():
                continue
            fields = line.split()
            if len(fields) < index:
                continue
            gene = fields[gene_col - 1]
            log2FoldChange = fields[lfc_col - 1]
            padj = fields[padj_col - 1]
            interactions.append((gene, log2FoldChange, padj))
    return interactions


print(load_deseq2_results(filename, gene_col, lfc_col, padj_col))
