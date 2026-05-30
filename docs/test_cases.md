# Casos de prueba

Estos casos permiten verificar si el programa funciona correctamente.

No son todavía pruebas automatizadas con pytest. Son escenarios para razonar, ejecutar y comparar resultados.

### Caso 1 · gen sobreexpresado significativo

Entrada:

```text
MX1	4.2	0.0001
```

Resultado esperado:

```text
El gen se reporta como significativo y sobreexpresado.
```

Criterio evaluado:

```text
El programa identifica correctamente genes con log2FoldChange positivo, padj significativo y magnitud suficiente.
```



### Caso 2 · gen subexpresado significativo

Entrada:

```text
GENE1	-3.0	0.001
```

Resultado esperado:

```text
El gen se reporta como significativo y subexpresado.
```

Criterio evaluado:

```text
El programa identifica correctamente genes con log2FoldChange negativo, padj significativo y magnitud suficiente.
```


### Caso 3 · gen no significativo por `padj`

Entrada:

```text
GENE2	3.0	0.8
```

Resultado esperado:

```text
El gen no aparece en el archivo de salida.
```

Criterio evaluado:

```text
El programa no debe reportar genes con padj alto, aunque tengan log2FoldChange grande.
```


### Caso 4 · gen no significativo por magnitud de cambio

Entrada:

```text
GENE3	0.3	0.001
```

Resultado esperado:

```text
El gen no aparece en el archivo de salida.
```

Criterio evaluado:

```text
El programa no debe reportar genes cuyo cambio sea pequeño, aunque tengan padj significativo.
```

## Caso: Sin archivo input

Entrada:

``` - o results/res.tsv
...
```

Resultado esperado:

```Error
...
```

Criterio evaluado: El archivo input no es ingresado


## Caso: Sin genes significativos

Entrada:

``` -i input.py - o results/res.tsv
...
```

Resultado esperado:

``` Resultados escritos en results/res.tsv
    Resumen del analisis:
    Total genes significativos: 0
    Upregulated genes: 0
    Downregulated genes: 0
...
```

Criterio evaluado: Filtro de significancia


## Caso: Sin genes significativos

Entrada:

``` -i input.py - o results/res.tsv
...
```

Resultado esperado:

``` Resultados escritos en results/res.tsv
    Resumen del analisis:
    Total genes significativos: 0
    Upregulated genes: 0
    Downregulated genes: 0
...
```

Criterio evaluado: Filtro de significancia

## Caso: Esperado

Entrada:

``` -i data/iav_deseq2_results.tsv -o results/res.tsv
...
```

Resultado esperado:

``` Resultados escritos en results/res.tsv
    Resumen del analisis:
    Total genes significativos: 100
    Upregulated genes: 60
    Downregulated genes: 40
...
```

Criterio evaluado: Funcionalidad del programa

### Caso límite 1 · valor `NA`

Entrada:

```text
GENE4	NA	0.001
```

Resultado esperado:

```text
La línea se ignora o se maneja con un mensaje claro.
El programa no debe romperse.
```

Criterio evaluado:

```text
El programa maneja valores no numéricos en columnas que deben convertirse a float.
```



### Caso límite 2 · línea incompleta

Entrada:

```text
GENE5	0.5
```

Resultado esperado:

```text
La línea se ignora.
El programa continúa procesando las demás líneas.
```

Criterio evaluado:

```text
El programa valida que existan suficientes columnas antes de intentar extraer datos.
```



### Caso límite 3 · archivo inexistente

Entrada:

```text
python analyze_iav.py data/no_existe.tsv results/iav_significant_genes.tsv
```

Resultado esperado:

```text
El programa muestra un mensaje claro indicando que el archivo no existe.
```

Criterio evaluado:

```text
El programa maneja errores de lectura del archivo de entrada.
```

## Tabla resumen de casos de prueba

| Caso | Entrada principal | Resultado esperado | Concepto evaluado |
|---|---|---|---|
| 1 | `MX1 4.2 0.0001` | aparece como `upregulated` | clasificación positiva |
| 2 | `GENE1 -3.0 0.001` | aparece como `downregulated` | clasificación negativa |
| 3 | `GENE2 3.0 0.8` | no aparece | filtro por `padj` |
| 4 | `GENE3 0.3 0.001` | no aparece | filtro por magnitud |
| 5 | `GENE4 NA 0.001` | se ignora | conversión numérica |
| 6 | `GENE5 0.5` | se ignora | línea incompleta |
| 7 | archivo inexistente | mensaje de error | manejo de errores |