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

    parser.add_argument("-g", "--gff-file", help="Ruta al archivo GFF de entrada")

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
    return sorted(interactions)


# ------------------------------------------
# Extraer la columna de atributos del archivo GFF
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: Extraer la columna de atributos del archivo GFF
# Entrada: gff_file
# Salida: gff_data (columna de atributos del archivo GFF)
# ------------------------------------------


def load_gff(gff_file):
    """Carga la columna de atributos de un archivo GFF.

    Lee un archivo GFF y extrae la columna de atributos (campo 9) de cada
    línea no comentada.

    Args:
        gff_file (str): Ruta al archivo GFF de entrada.

    Returns:
        list[str]: Lista de cadenas con la columna de atributos de cada línea.

    Raises:
        SystemExit: Si el archivo no existe.
    """

    gff_data = []
    if not os.path.exists(gff_file):
        print("Error: archivo no encontrado")
        exit(1)
    index = 9
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < index:
                continue
            attributes = fields[8]
            gff_data.append(attributes)
    return gff_data


# ------------------------------------------
# Extraer los datos significativos del archivo GFF
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: Extraer los datos significativos del archivo GFF
# Entrada: gff_data
# Salida: Informacion significativa del archivo GFF ordenada
# ------------------------------------------


def save_gff_data(gff_data):
    """Extrae campos significativos de la columna de atributos del GFF.

    Busca atributos `ID`, `Name`, `description` y `gene_type` y devuelve una
    lista de tuplas ordenadas con (gene_id, gene_name, description, gene_type).

    Args:
        gff_data (iterable): Iterable de cadenas (atributos) extraídas del GFF.

    Returns:
        list[tuple[str, str, str, str]]: Lista ordenada de tuplas con los campos.
    """

    all_data = []
    for attributes in gff_data:
        if (
            "ID" in attributes
            and "Name" in attributes
            and "description" in attributes
            and "gene_type" in attributes
        ):
            gene_id = attributes.split("ID=")[1].split(";")[0]
            gene_name = attributes.split("Name=")[1].split(";")[0]
            description = attributes.split("description=")[1].split(";")[0]
            gene_type = attributes.split("gene_type=")[1].split(";")[0]
            all_data.append(
                (
                    gene_id,
                    gene_name,
                    description,
                    gene_type,
                )
            )
    return sorted(all_data)


# ------------------------------------------
# Combinar los datos del archivo de entrada con los datos del archivo GFF
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: Combinar los datos del archivo de entrada con los datos del archivo GFF
# Entrada: genes_tupla, gff_significant_data
# Salida: merged_data (datos combinados)
# ------------------------------------------


def merge_data(genes_tupla, gff_significant_data):
    """Combina resultados de DESeq2 con datos significativos del GFF.

    Busca coincidencias entre el nombre del gen en `genes_tupla` y `gene_name`
    en `gff_significant_data` y produce una lista de tuplas combinadas.

    Args:
        genes_tupla (iterable): Iterable de tuplas (gene, log2FoldChange, padj).
        gff_significant_data (iterable): Iterable de tuplas
            (gene_id, gene_name, description, gene_type).

    Returns:
        list[tuple]: Lista de tuplas combinadas
            (gene, log2FoldChange, padj, gene_id, gene_name, description, gene_type).
    """

    merged = []
    for gene, log2FoldChange, padj in genes_tupla:
        for gene_id, gene_name, description, gene_type in gff_significant_data:
            if gene == gene_name:
                merged.append(
                    (
                        gene,
                        log2FoldChange,
                        padj,
                        gene_id,
                        gene_name,
                        description,
                        gene_type,
                    )
                )
                break
    return merged


# ------------------------------------------
# Extraer los genes significativos del archivo de entrada
# ------------------------------------------
# ------------------------------------------
# Responsabilidad: Guardar los genes significativos del archivo de entrada
# Entrada: genes_tupla
# Salida: filtro (genes significativos)
# ------------------------------------------


def is_significant(merged_data, padj_threshold, lfc_threshold):
    """Filtra genes significativos según umbrales de padj y log2FoldChange.

    Args:
        merged_data (iterable): Iterador o lista de tuplas `(gene, log2FoldChange, padj, gene_id, gene_name, description, gene_type)`.
        padj_threshold (float): Umbral máximo para `padj` (valor < padj_threshold).
        lfc_threshold (float): Umbral mínimo absoluto para `log2FoldChange`.

    Returns:
        list[tuple[str, str, str, str, str, str, str]]: Subconjunto de `merged_data` que cumplen los umbrales.
    """

    filtro = []
    for (
        gene,
        log2FoldChange,
        padj,
        gene_id,
        gene_name,
        description,
        gene_type,
    ) in merged_data:
        if abs(float(log2FoldChange)) >= lfc_threshold and float(padj) < padj_threshold:
            filtro.append(
                (gene, log2FoldChange, padj, gene_id, gene_name, description, gene_type)
            )
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
    for (
        gene,
        log2FoldChange,
        padj,
        gene_id,
        gene_name,
        description,
        gene_type,
    ) in filtro:
        if float(log2FoldChange) >= 1:
            classify.append(
                (
                    gene,
                    log2FoldChange,
                    padj,
                    gene_id,
                    gene_name,
                    description,
                    gene_type,
                    "upregulated",
                )
            )
        elif float(log2FoldChange) <= -1:
            classify.append(
                (
                    gene,
                    log2FoldChange,
                    padj,
                    gene_id,
                    gene_name,
                    description,
                    gene_type,
                    "downregulated",
                )
            )
        else:
            classify.append(
                (
                    gene,
                    log2FoldChange,
                    padj,
                    gene_id,
                    gene_name,
                    description,
                    gene_type,
                    "not_significant",
                )
            )
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
        classify (iterable): Lista de tuplas `(gene, log2FoldChange, padj, gene_id, gene_name, description, gene_type, classification)`.

    Returns:
        None
    """

    with open(output, "w") as f:
        f.write(
            "Gene\tlog2FoldChange\tpadj\tGene_ID\tGene_Name\tDescription\tGene_Type\tClassification\n"
        )
        for (
            gene,
            log2FoldChange,
            padj,
            gene_id,
            gene_name,
            description,
            gene_type,
            classification,
        ) in classify:
            f.write(
                f"{gene}\t{log2FoldChange}\t{padj}\t{gene_id}\t{gene_name}\t{description}\t{gene_type}\t{classification}\n"
            )
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
        1
        for _, _, _, _, _, _, _, classification in classify
        if classification == "upregulated"
    )
    downregulated = sum(
        1
        for _, _, _, _, _, _, _, classification in classify
        if classification == "downregulated"
    )
    not_significant = sum(
        1
        for _, _, _, _, _, _, _, classification in classify
        if classification == "not_significant"
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
    gff_file = args.gff_file
    padj_threshold = args.padj_threshold
    lfc_threshold = args.lfc_threshold
    gene_col = args.columna_gene
    lfc_col = args.columna_lfc
    padj_col = args.columna_padj

    genes_tupla = load_deseq2_results(filename, gene_col, lfc_col, padj_col)

    gff_data = load_gff(gff_file)

    gff_significant_data = save_gff_data(gff_data)

    merged_data = merge_data(genes_tupla, gff_significant_data)

    filtro = is_significant(merged_data, padj_threshold, lfc_threshold)

    classify = sorted(classify_gene(filtro))

    write_results(output, classify)

    print_summary(classify)


if __name__ == "__main__":
    main()
