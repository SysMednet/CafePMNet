# CafePMNets

CafePMNets is a multi-scale network that modelling protein-protein, module-module and protein-module interactions. The construction of CafePMNets is comprised of three main steps:  
* Step 1: Protein-protein interaction network construction 
* Step 2: Module-module interaction network construction 
* Step 3: Protein-module interaction network construction

# Input File Format

1.	Experimental protein-protein interactions file (tab-delimited; with header)
    - Column 1: Protein 1; Column 2: Protein 2
2.	Co-expressed gene pairs file (tab-delimited; no header):
    - Column 1: Protein 1; Column 2: Protein 2; Column 3: Pearson correlation coefficient
3.	BP semantic similarity score file (tab-delimited; no header):
    - Column 1: Protein 1; Column 2: Protein 2; Column 3: BP score
4.	CC semantic similarity score file (tab-delimited; no header):
    - Column 1: Protein 1; Column 2: Protein 2; Column 3: CC score
5.	Reference genome file (tab-delimited; no header):
    - Column 1: Entry (Accession); Column 2: Entry name (ID); Column 3: Gene names
    - Column 4: Protein names; Column 5: Organism ID; Column 6: Status
6.	Modules file (tab-delimited; with header):
    - Column 1: Module ID; Column 2: Number of module protein; Column 3: Module protein (Uniprot AC)

# Basic Usage

Run the main script with your expression and weight files. Note that `utils.py` must be in the same directory as the main script to ensure internal functions load correctly.

**Step 1: Protein-protein interaction network (PPIN) construction**

The example dataset is stored in the Dataset folder, all of the example outputs are provided in the PPIN folder.

Step 1-1: Calculate the $R^2$ of PPIN under different PCC thresholds.

```bash
python3 1-1_R2_PPIN.py –e ./Dataset/9606_expPPI_network_OnlyAC.txt –p ./Dataset/COAD_PCC.txt –o ./PPIN/R2_COAD_PPIN.txt
```
`-h`: Get help with the commands.  
`-e`: Experimental PPI file.  
`-p`: Gene co-expression file.  
`-o`: Output file containing average $R^2$ values and edge counts for each PCC threshold.

Step 1-2: Calculate the average RSS of PPIN under different PCC thresholds.

```bash
python3 1-2_AvgRSS_PPIN.py –e ./Dataset/9606_expPPI_network_ OnlyAC.txt –p ./Dataset/COAD_PCC.txt –b ./Dataset/ALL_PPI_BPscore.txt –c ./Dataset/ALL_PPI_CCscore.txt –o ./PPIN/AvgRSS_COAD_PPIN.txt
```
`-h`: Get help with the commands.  
`-e`: Experimental PPI file.  
`-p`: Gene co-expression file.  
`-b`: BP semantic similarity score file.  
`-c`: CC semantic similarity score file.  
`-o`: Output file containing average RSS values and edge counts for each PCC threshold.  
`--threshold_start`: Starting PCC threshold (default: 0.1).  
`--threshold_end`: Ending PCC threshold (default: 0.9).  
`--threshold_step`: PCC threshold increment (default: 0.1).

Step 1-3: Calculate the geometric mean scores and plotting of PPIN under different PCC thresholds.

```bash
python3 1-3_geometric_mean_plotting.py –r ./PPIN/R2_COAD_PPIN.txt –s ./PPIN/AvgRSS_COAD_PPIN.txt –o ./PPIN/GM_COAD_PPIN.txt
```
`-h`: Get help with the commands.  
`-r`: R² result file.  
`-s`: Average RSS result file.  
`-o`: Output file containing geometric mean scores. A threshold selection plot with the same basename will also be generated.

Step 1-4: Construct PPIN with assigned threshold.

```bash
python3 1-4_construct_PPIN.py –e ./Dataset/9606_expPPI_network_ OnlyAC.txt –p ./Dataset/COAD_PCC.txt –t 0.5 –o ./PPIN/COAD_PPIN.txt
```
`-h`: Get help with the commands.  
`-e`: Experimental PPI file.  
`-p`: Gene co-expression file.  
`-t`: PCC threshold for PPIN construction.  
`-o`: Output PPIN file.

**Step 2: Module-module interaction network (MMIN) construction**

The required datasets are stored in the Dataset folder, all of the example outputs are provided in the MMIN folder.

Step 2-1: Calculate hypergeometric p-value based on experimental PPI or co-expressed gene pairs for each module-module pair.

```bash
python3 2-1_cor_expPPI_MMI_hyperP.py --mode coexpressed –m ./Dataset/ Modules.txt –n ./Dataset/COAD_PCC.txt --threshold 0.5 –o ./MMIN/coexp_MMIN_hyperP.txt
```
```bash
python3 2-1_cor_expPPI_ MMI_hyperP.py --mode expPPI –m ./Dataset/Modules.txt –n ./Dataset/9606_expPPI_network_ OnlyAC.txt --threshold 0.5 –o ./MMIN/expPPI_MMIN_hyperP.txt
```
`-help`: Get help with the commands.  
`-m`: Module file.  
`-n`: Interaction network file.  
`-o`: Output file containing hypergeometric enrichment results.  
`--mode`: Reference interaction source ("expPPI" or "coexpressed").  
`--threshold`: Absolute PCC threshold (required for coexpressed mode).

Step 2-2: Calculate the $R^2$ of MMIN under different meta-z thresholds.

```bash
python3 2-2_R2_MMIN.py –e ./MMIN/expPPI_MMIN_hyperP.txt –p ./MMIN/coexp_MMIN_hyperP.txt –o ./MMIN/R2_COAD_MMIN.txt
```
`-h`: Get help with the commands.  
`-e`: MMI hypergeometric enrichment result file based on experimental PPI.  
`-p`: MMI hypergeometric enrichment result file based on gene co-expression.  
`-o`: Output file for the calculated $R^2$ values across all meta-z thresholds (can be .txt or .csv).  
`--threshold_start`: Starting meta-z threshold (default: 1.0).  
`--threshold_end`: Ending meta-z threshold (default: 9.0).  
`--threshold_step`: Meta-z threshold increment (default: 0.5).

Step 2-3: Calculate the joint average RSS of MMIN under different meta-z thresholds.

```bash
python3 2-3_joint_AvgRSS_MMIN.txt –k ./Dataset/9606_expPPI_network_OnlyAC.txt –e ./MMIN/expPPI_MMIN_hyperP.txt –b ./Dataset/ALL_PPI_BPscore.txt –c ./Dataset/ALL_PPI_CCscore.txt –m ./Dataset/Modules.txt –p ./MMIN/coexp_MMIN_hyperP.txt –o ./MMIN/joint_AvgRSS_COAD_MMIN.txt
```
`-k`: Experimental PPI file.  
`-e`: MMI hypergeometric enrichment result file based on experimental PPI.  
`-b`: BP semantic similarity score file.  
`-c`: CC semantic similarity score file.  
`-m`: Module file.  
`-p`: MMI hypergeometric enrichment result file based on co-expression.  
`-o`: Output file containing joint average RSS values for each meta-z threshold (can be .txt or .csv).  
`--threshold_start`: Starting meta-z threshold (default: 1.0).  
`--threshold_end`: Ending meta-z threshold (default: 9.0).  
`--threshold_step`: Meta-z threshold increment (default: 0.5).

Step 2-4: Calculate the geometric mean scores and plotting of MMIN under different meta-z thresholds.

```bash
python3 2-4_geometric_mean_plotting.py –r ./MMIN/R2_COAD_MMIN.txt –s ./MMIN/joint_AvgRSS_COAD_MMIN.txt –o ./MMIN/GM_COAD_MMIN.txt
```
`-h`: Get help with the commands.  
`-r`: R² result file.  
`-s`: Average RSS result file.  
`-o`: Output file containing geometric mean scores. A threshold selection plot with the same basename will also be generated.

Step 2-5: Construct MMIN with assigned meta-z threshold.

```bash
python3 2-5_construct_MMIN.py –e ./MMIN/expPPI_MMIN_hyperP.txt –p ./MMIN/coexp_MMIN_hyperP.txt –t 7.5 –o ./MMIN/COAD_MMIN.txt
```
`-h`: Get help with the commands.  
`-e`: MMI hypergeometric enrichment result file based on experimental PPI.  
`-p`: MMI hypergeometric enrichment result file based on gene co-expression.  
`-t`: Meta-z threshold for MMIN construction.  
`-o`: Output MMIN file.

**Step 3: Protein-module interaction network (CafePMNets) construction**

Step 3-1: Calculate hypergeometric p-value based on experimental PPI or co-expressed gene pairs for each protein-module pair.

```bash
python3 3-1_cor_ expPPI_PMI_hyperP.py --mode coexpressed –m ./Dataset/Modules.txt –n ./Dataset/COAD_PCC.txt –k ./Dataset/9606_Uniprot_Reference_Genome.txt --threshold 0.5 –o ./PMIN/coexp_PMI_hyperP.txt
```
```bash
python3 3-1_cor_expPPI_ PMI_hyperP.py --mode expPPI –m ./Dataset/Modules.txt –n ./Dataset/9606_expPPI_network_OnlyAC.txt –k ./Dataset/9606_Uniprot_Reference_Genome.txt --threshold 0.5 –o ./PMIN/expPPI_PMI_hyperP.txt
```

`-help`: Get help with the commands.  
`-m`: Module file.  
`-n`: Interaction network file.  
`-k`: Reference genome file.  
`-o`: Output file containing hypergeometric enrichment results.  
`--mode`: Reference interaction source ("expPPI" or "coexpressed").  
`--threshold`: Absolute PCC threshold (required for coexpressed mode).

Step 3-2: Calculate the $R^2$ of CafePMNets under different meta-z thresholds.

```bash
python3 3-2_R2_CafePMNets.py –e ./PMIN/expPPI_PMIN_hyperP.txt –p ./PMIN/coexp_PMIN_hyperP.txt –pp ./PPIN/COAD_PPIN.txt –mm ./MMIN/COAD_MMIN.txt -o ./PMIN/R2_COAD_CafePMNets.txt
```

`-h`: Get help with the commands.  
`-e`: PMI hypergeometric enrichment result file based on experimental PPI.  
`-p`: PMI hypergeometric enrichment result file based on gene co-expression.  
`-pp`: PPIN file.  
`-mm`: MMIN file.  
`-o`: Output file for the calculated $R^2$ values across all meta-z thresholds (can be .txt or .csv).  
`--threshold_start`: Starting meta-z threshold (default: 1.0).  
`--threshold_end`: Ending meta-z threshold (default: 9.0).  
`--threshold_step`: Meta-z threshold increment (default: 0.5).

Step 3-3: Calculate the joint average RSS of CafePMNets under different meta-z thresholds.

```bash
python3 3-3_joint_AvgRSS_CafePMNets.txt –e ./PMIN/expPPI_PMIN_hyperP.txt –b ./Dataset/ALL_PPI_BPscore.txt –c ./Dataset/ALL_PPI_CCscore.txt –m ./Dataset/Modules.txt –p ./PMIN/coexp_PMIN_hyperP.txt –pp ./PPIN/COAD_PPIN.txt –mm ./MMIN/COAD_MMIN.txt –o ./PMIN/joint_AvgRSS_COAD_CafePMNets.txt
```

`-e`: PMI hypergeometric enrichment result file based on experimental PPI.  
`-b`: BP semantic similarity score file.  
`-c`: CC semantic similarity score file.  
`-m`: Module file.  
`-p`: PMI hypergeometric enrichment result file based on co-expression.  
`-pp`: PPIN file.  
`-mm`: MMIN file.  
`-o`: Output file containing joint average RSS values for each meta-z threshold (can be .txt or .csv).  
`--threshold_start`: Starting meta-z threshold (default: 1.0).  
`--threshold_end`: Ending meta-z threshold (default: 9.0).  
`--threshold_step`: Meta-z threshold increment (default: 0.5).

Step 3-4: Calculate the geometric mean scores and plotting of CafePMNets under different meta-z thresholds.

```bash
python3 3-4_geometric_mean_plotting.py –r ./PMIN/R2_COAD_CafePMNets.txt –s ./PMIN/joint_AvgRSS_COAD_CafePMNets.txt –o ./PMIN/GM_COAD_CafePMNets.txt
```

`-h`: Get help with the commands.  
`-r`: R² result file.  
`-s`: Average RSS result file.  
`-o`: Output file containing geometric mean scores. A threshold selection plot with the same basename will also be generated.

Step 3-5: Construct CafePMNets with assigned meta-z threshold.

```bash
python3 3-5_construct_CafePMNets.py –e ./PMIN/expPPI_PMI_hyperP.txt –p ./PMIN/coexp_PMI_hyperP.txt –pp ./PPIN/COAD_PPIN.txt –mm ./MMIN/COAD_MMIN.txt –t 4.5 –o ./PMIN/COAD_CafePMNets.txt
```

`-h`: Get help with the commands.  
`-e`: PMI hypergeometric enrichment result file based on experimental PPI.  
`-p`: PMI hypergeometric enrichment result file based on gene co-expression.  
`-pp`: PPIN file.  
`-mm`: MMIN file.  
`-t`: Meta-z threshold for CafePMNets construction.  
`-o`: Output CafePMNets file.

# Dependencies

The code was developed and tested with the following software:

- Python 3.12.4
- NumPy 1.26.4
- pandas 2.2.2
- SciPy 1.13.1
- NetworkX 3.2.1
- matplotlib 3.8.4
- seaborn 0.13.2

Install all dependencies using:

```bash
pip install -r requirements.txt
```
