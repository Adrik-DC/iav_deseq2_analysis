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
            `input`, `output`, `columna_gene`, `columna_lfc`, `columna_padj`.
    """

    parser = argparse.ArgumentParser(description=" Analisis de IAV")

    parser.add_argument("-i", "--input", help="Ruta al archivo fasta de entrada")

    parser.add_argument("-o", "--output", help="Ruta al archivo fasta de salida")

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
output = argumentos.output
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


# print(load_deseq2_results(filename, gene_col, lfc_col, padj_col))

genes_tupla = load_deseq2_results(filename, gene_col, lfc_col, padj_col)

# ------------------------------------------
# Extraer los genes significativos del archivo de entrada
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: Guardar los genes significativos del archivo de entrada
# Entrada: genes_tupla
# Salida: filtro (genes significativos)
# ------------------------------------------


def is_significant(genes_tupla):
    filtro = []
    for gene, log2FoldChange, padj in genes_tupla:
        if abs(float(log2FoldChange)) >= 1 and float(padj) < 0.05:
            filtro.append((gene, log2FoldChange, padj))
    return filtro


filtro = is_significant(genes_tupla)
# print(filtro)

# ------------------------------------------
# Filtrar los genes por regulación
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: filtrar los genes por regulación del archivo de entrada
# Entrada: filttro (genes significativos)
# Salida: classify (genes regulados)
# ------------------------------------------


def classify_gene(filtro):
    classify = []
    for gene, log2FoldChange, padj in filtro:
        if float(log2FoldChange) >= 1:
            classify.append((gene, log2FoldChange, padj, "upregulated"))
        elif float(log2FoldChange) <= -1:
            classify.append((gene, log2FoldChange, padj, "downregulated"))
    return classify


classify = classify_gene(filtro)
# print(classify)

# ------------------------------------------
# Escribir los resultados en un archivo de salida (.tsv)
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: Escribir los resultados en un archivo de salida (.tsv)
# Entrada: output, classify
# Salida: output (archivo de salida con los resultados)
# ------------------------------------------


def write_results(output, classify):
    with open(output, "w") as f:
        f.write("Gene\tlog2FoldChange\tpadj\tClassification\n")
        for gene, log2FoldChange, padj, classification in classify:
            f.write(f"{gene}\t{log2FoldChange}\t{padj}\t{classification}\n")
    print(f"Resultados escritos en {output}")


write_results(output, classify)

# ------------------------------------------
# Escribir el resumen del análisis en la consola
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: Escribir el resumen del análisis en la consola
# Entrada: classify
# Salida: None
# ------------------------------------------


def print_summary(classify):
    upregulated = sum(
        1 for _, _, _, classification in classify if classification == "upregulated"
    )
    downregulated = sum(
        1 for _, _, _, classification in classify if classification == "downregulated"
    )
    print("Resumen del análisis:")
    print(f"Total genes significativos: {len(classify)}")
    print(f"Upregulated genes: {upregulated}")
    print(f"Downregulated genes: {downregulated}")


print_summary(classify)
