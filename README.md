# IAV DESeq2 Analysis

Programa en Python para identificar genes diferencialmente expresados a partir de resultados de DESeq2 durante infección por Influenza A Virus.
```



### 2. Estructura del proyecto

Ejemplo:

```text
iav_deseq2_analysis/
├── analyse_iav.py
├── data
│   ├── human_genes.gff
│   └── iav_deseq2_results.tsv
├── docs
│   ├── context.md
│   ├── design.md
│   └── test_cases.md
├── pyproject.toml
├── README.md
├── results
│   └── res.tsv
```


### 3. Cómo ejecutar el programa

Ejemplo:

```bash
uv run python3 analyse_iav.py -i data/iav_deseq2_results.tsv -o results/res.tsv

```


### 4. Uso de thresholds opcionales

Ejemplo:

```bash
uv run python3 analyse_iav.py -i data/iav_deseq2_results.tsv -o results/res.tsv --lfc-threshold 0.0 --padj-threshold -1 --columna-gene 1 --columna-lfc 3 --columna-padj 7
```



### 5. Valores por defecto

Documentar que:

- `--lfc_threshold` usa `1.0`
- `--padj_threshold` usa `0.05`
- `--columna-gene` usa 1
- `--columna-lfc` usa 3
- `--columna-padj` usa 7

si el usuario no proporciona argumentos.



### 6. Descripción breve de salida esperada

Ejemplo:

```text
gene	log2FoldChange	padj	status
MX1	4.2	0.0001	upregulated
```

## Actualizacion para agregar informacion del archivo gff

El programa ahora agrega informacion del gen del archivo .gff

- `--gff-file`

### Cómo ejecutar el programa (actualizacion)

Ejemplo:

```bash
uv run python3 analyse_iav.py -i data/iav_deseq2_results.tsv -o results/res.tsv -g data/human_genes.gff 
```

### Uso de thresholds opcionales (actualizacion)

Ejemplo:

```bash
uv run python3 analyse_iav.py -i data/iav_deseq2_results.tsv -o results/res.tsv -g data/human_genes.gff --lfc-threshold 0.0 --padj-threshold -1 --columna-gene 1 --columna-lfc 3 --columna-padj 7
```


## Ejemplo de salida esperada (actualizacion)

```text
gene	log2FoldChange	padj    Gene_ID	Gene_Name   Description Gene_Type   Classification
MX1	4.2	0.0001  ENSG00000278540_MX1	 MX1  acetyl-CoA carboxylase alpha, first step...	upregulated
IFIT1	5.1	0.00001  ENSG00000131473_IFIT1	IFIT1	ATP-citrate lyase, links TCA cycle...  upregulated
```