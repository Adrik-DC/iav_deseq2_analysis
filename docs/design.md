# Diseño del programa

## Idea general

El programa leerá un archivo de resultados de DESeq2 línea por línea. Para cada gen, extraerá el nombre, log2FoldChange y padj. Después evaluará si el gen cumple los criterios de significancia y lo clasificará como sobreexpresado o subexpresado.

## Estructuras de datos propuestas

Se usará una lista para almacenar genes significativos.

Cada gen leído puede representarse como una tupla:

(gene, log2FoldChange, padj)

Cada gen significativo puede representarse como una tupla:

(gene, log2FoldChange, padj, status)

Donde status puede ser:

- upregulated
- downregulated

## Algoritmo general

1. Leer argumentos de entrada y salida.
2. Abrir el archivo de entrada.
3. Crear una lista vacía para guardar genes leídos.
4. Leer el archivo línea por línea.
5. Ignorar encabezado, líneas vacías o líneas incompletas.
6. Extraer gene, log2FoldChange y padj.
7. Convertir log2FoldChange y padj a números.
8. Guardar cada gen válido en una lista.
9. Crear una lista vacía para genes significativos.
10. Recorrer los genes válidos.
11. Evaluar si cada gen cumple padj < 0.05 y abs(log2FoldChange) >= 1.
12. Si cumple, clasificarlo como sobreexpresado o subexpresado.
13. Guardarlo en la lista de resultados.
14. Escribir la lista en el archivo de salida.
15. Imprimir un resumen final.

## Funciones sugeridas
1. parsear_argumentos()
   - Lee los arguentos mediante argparse

2. load_deseq2_results(filename)
   - Lee el archivo de entrada y regresa una lista de genes válidos.

3. is_significant(genes_tupla)
   - Evalúa si un gen cumple los criterios de significancia.

4. classify_gene(filtro)
   - Clasifica el gen como sobreexpresado o subexpresado.

5. write_results(output_file, classify)
   - Guarda los resultados filtrados.

6. print_summary(classify)
   - Muestra el resumen final.

7. main()
   - Coordina el flujo general del programa.

flowchart TD
    A[Inicio] --> B[Leer argumentos]
    B --> C[Abrir archivo TSV]
    C --> D[Leer líneas del archivo]
    D --> E{¿La línea es válida?}
    E -->|No| D
    E -->|Sí| F[Extraer gene, log2FC y padj]
    F --> G[Convertir log2FC y padj a números]
    G --> H{¿padj < 0.05 y abs(log2FC) >= 1?}
    H -->|No| D
    H -->|Sí| I{¿log2FC > 0?}
    I -->|Sí| J[Guardar como upregulated]
    I -->|No| K[Guardar como downregulated]
    J --> D
    K --> D
    D --> L{¿Terminó el archivo?}
    L -->|No| D
    L -->|Sí| M[Escribir resultados]
    M --> N[Mostrar resumen]
    N --> O[Fin]

    