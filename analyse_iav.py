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
            `input`, `output`, `columna_gene`, `columna_lfc`, `columna_padj`,
            `lfc_threshold`, `padj_threshold`.
    """

    parser = argparse.ArgumentParser(description=" Analisis de IAV")

    parser.add_argument("-i", "--input", help="Ruta al archivo fasta de entrada")

    parser.add_argument("-o", "--output", help="Ruta al archivo fasta de salida")

    parser.add_argument(
        "--lfc-threshold",
        type=float,
        default=1,
        help="Umbral para el log2FoldChange",
    )

    parser.add_argument(
        "--padj-threshold",
        type=float,
        default=0.05,
        help="Umbral para el padj (>=0 && <1)",
    )

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


def load_deseq2_results(filename, gene_col, lfc_col, padj_col):
    """Carga los resultados de DESeq2 desde un archivo de texto.

    Lee el archivo de resultados (esperando una cabecera en la primera línea)
    y extrae las columnas de interés para cada gen.

    Args:
        filename (str): Ruta al archivo de entrada.
        gene_col (int): Índice de la columna que contiene el nombre del gen.
        lfc_col (int): Índice de la columna que contiene el log2FoldChange.
        padj_col (int): Índice de la columna que contiene el padj.

    Returns:
        list[tuple[str, str, str]]: Lista de tuplas `(gene, log2FoldChange, padj)`.

    Raises:
        SystemExit: Si el archivo no existe.
    """

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


# ------------------------------------------
# Extraer los genes significativos del archivo de entrada
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: Guardar los genes significativos del archivo de entrada
# Entrada: genes_tupla
# Salida: filtro (genes significativos)
# ------------------------------------------


def is_significant(genes_tupla, padj_threshold, lfc_threshold):
    """Filtra genes significativos según umbrales de padj y log2FoldChange.

    Args:
        genes_tupla (iterable): Iterador o lista de tuplas `(gene, log2FoldChange, padj)`.
        padj_threshold (float): Umbral máximo para `padj` (valor < padj_threshold).
        lfc_threshold (float): Umbral mínimo absoluto para `log2FoldChange`.

    Returns:
        list[tuple[str, str, str]]: Subconjunto de `genes_tupla` que cumplen los umbrales.
    """

    filtro = []
    for gene, log2FoldChange, padj in genes_tupla:
        if abs(float(log2FoldChange)) >= lfc_threshold and float(padj) < padj_threshold:
            filtro.append((gene, log2FoldChange, padj))
    return filtro


# ------------------------------------------
# Filtrar los genes por regulación
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: filtrar los genes por regulación del archivo de entrada
# Entrada: filttro (genes significativos)
# Salida: classify (genes regulados)
# ------------------------------------------


def classify_gene(filtro):
    """Clasifica genes según la dirección de regulación.

    Args:
        filtro (iterable): Lista o iterador de tuplas `(gene, log2FoldChange, padj)`
            que ya han sido filtradas como significativas.

    Returns:
        list[tuple[str, str, str, str]]: Lista de tuplas
            `(gene, log2FoldChange, padj, classification)` donde
            `classification` es una de `"upregulated"`, `"downregulated"`,
            o `"not_significant"`.
    """

    classify = []
    for gene, log2FoldChange, padj in filtro:
        if float(log2FoldChange) >= 1:
            classify.append((gene, log2FoldChange, padj, "upregulated"))
        elif float(log2FoldChange) <= -1:
            classify.append((gene, log2FoldChange, padj, "downregulated"))
        else:
            classify.append((gene, log2FoldChange, padj, "not_significant"))
    return classify


# ------------------------------------------
# Escribir los resultados en un archivo de salida (.tsv)
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: Escribir los resultados en un archivo de salida (.tsv)
# Entrada: output, classify
# Salida: output (archivo de salida con los resultados)
# ------------------------------------------


def write_results(output, classify):
    """Escribe los resultados clasificados en un archivo TSV.

    Args:
        output (str): Ruta del archivo de salida donde se escribirán los resultados.
        classify (iterable): Lista de tuplas `(gene, log2FoldChange, padj, classification)`.

    Returns:
        None
    """

    with open(output, "w") as f:
        f.write("Gene\tlog2FoldChange\tpadj\tClassification\n")
        for gene, log2FoldChange, padj, classification in classify:
            f.write(f"{gene}\t{log2FoldChange}\t{padj}\t{classification}\n")
    print(f"Resultados escritos en {output}")


# ------------------------------------------
# Escribir el resumen del análisis en la consola
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: Escribir el resumen del análisis en la consola
# Entrada: classify
# Salida: None
# ------------------------------------------


def print_summary(classify):
    """Imprime un resumen estadístico simple del análisis en consola.

    Args:
        classify (iterable): Lista de tuplas `(gene, log2FoldChange, padj, classification)`.

    Returns:
        None
    """

    upregulated = sum(
        1 for _, _, _, classification in classify if classification == "upregulated"
    )
    downregulated = sum(
        1 for _, _, _, classification in classify if classification == "downregulated"
    )
    not_significant = sum(
        1 for _, _, _, classification in classify if classification == "not_significant"
    )
    print("Resumen del análisis:")
    print(f"Total genes significativos: {len(classify)}")
    print(f"Upregulated genes: {upregulated}")
    print(f"Downregulated genes: {downregulated}")
    print(f"Genes no significativos: {not_significant}")


def main():
    """Punto de entrada principal del script.

    Ejecuta el flujo completo: parsear argumentos, cargar datos, filtrar,
    clasificar, escribir resultados e imprimir resumen.

    Returns:
        None
    """

    args = parsear_argumentos()
    filename = args.input
    output = args.output
    padj_threshold = args.padj_threshold
    lfc_threshold = args.lfc_threshold
    gene_col = args.columna_gene
    lfc_col = args.columna_lfc
    padj_col = args.columna_padj

    genes_tupla = load_deseq2_results(filename, gene_col, lfc_col, padj_col)

    filtro = is_significant(genes_tupla, padj_threshold, lfc_threshold)

    classify = classify_gene(filtro)

    write_results(output, classify)

    print_summary(classify)


if __name__ == "__main__":
    main()
