# Casos de prueba

Estos casos permiten verificar si el programa funciona correctamente.

No son todavía pruebas automatizadas con pytest. Son escenarios para razonar, ejecutar y comparar resultados.

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
