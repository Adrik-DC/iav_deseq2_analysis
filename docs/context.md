# Contexto del proyecto

## Problema

Se desea analizar un archivo de resultados de DESeq2 para identificar genes diferencialmente expresados durante una infección por Influenza A Virus.

## Objetivo del programa

Construir una primera versión de un programa en Python que lea un archivo TSV, filtre genes significativos y clasifique genes sobreexpresados y subexpresados.

## Archivos de entrada

- data/iav_deseq2_results.tsv

## Archivo de salida esperado

- results/iav_significant_genes.tsv

## Ejemplo de salida esperada

```text
gene	log2FoldChange	padj	status
MX1	4.2	0.0001	upregulated
IFIT1	5.1	0.00001	upregulated
```

## Requisitos funcionales

1. Leer un archivo TSV con resultados de DESeq2.
2. Extraer el nombre del gen, log2FoldChange y padj.
3. Ignorar líneas inválidas o incompletas.
4. Identificar genes significativos usando padj < 0.05 y abs(log2FoldChange) >= 1.
5. Clasificar genes como sobreexpresados o subexpresados.
6. Mostrar un resumen en pantalla.
7. Guardar los genes filtrados en un archivo de salida.

## Criterios de significancia

Un gen se considerará diferencialmente expresado si cumple:

- padj < 0.05
- abs(log2FoldChange) >= 1

Si log2FoldChange > 0, el gen se clasificará como sobreexpresado.
Si log2FoldChange < 0, el gen se clasificará como subexpresado.

## Nueva Funcionalidad

El programa ahora permite configurar los thresholds desde línea de comandos mediante:

- `--lfc_threshold`
- `--padj_threshold`

Si el usuario no proporciona estos argumentos, el programa utiliza:

- `1.0` para log2 Fold Change
- `0.05` para adjusted p-value

## Actualizacion para agregar informacion del archivo gff

El programa ahora agrega informacion del gen del archivo .gff

- `--gff-file`

## Archivos de entrada (actualizacion)

- data/iav_deseq2_results.tsv

## Archivo de salida esperado (actualizacion)

- results/iav_significant_genes.tsv

## Archivo de gff esperado (actualizacion)

- data/human_genes.gff

## Ejemplo de salida esperada (actualizacion)

```text
gene	log2FoldChange	padj    Gene_ID	Gene_Name   Description Gene_Type   Classification
MX1	4.2	0.0001  ENSG00000278540_MX1	 MX1  acetyl-CoA carboxylase alpha, first step...	upregulated
IFIT1	5.1	0.00001  ENSG00000131473_IFIT1	IFIT1	ATP-citrate lyase, links TCA cycle...  upregulated
```

## Requisitos funcionales (actualizacion)

1. Leer un archivo TSV con resultados de DESeq2 y archivo gff con informacion de los genes.
2. Extraer el nombre del gen, log2FoldChange, padj, gene_id, gene_name, description, gene_type.
3. Ignorar líneas inválidas o incompletas.
4. Identificar genes significativos usando padj < 0.05 y abs(log2FoldChange) >= 1.
5. Clasificar genes como sobreexpresados o subexpresados.
6. Mostrar un resumen en pantalla.
7. Guardar los genes filtrados en un archivo de salida.