# CafePMNets

CafePMNets is a multi-scale network that modelling protein-protein, module-module and protein-module interactions. The construction of CafePMNets is comprised of three main steps:  
* Step 1: Protein-protein interaction network construction 
* Step 2: Module-module interaction network construction 
* Step 3: Protein-module interaction network construction

# Input File Format

1.	Experimental protein-protein interactions file (tab-delimited; with header)
  - Column 1: Protein 1
  - Column 2: Protein 2
2.	Co-expressed gene pairs file (tab-delimited; no header):
  - Column 1: Protein 1
  - Column 2: Protein 2
  - Column 3: Pearson correlation coefficient
3.	BP semantic similarity score file (tab-delimited; no header):
  - Column 1: Protein 1
  - Column 2: Protein 2
  - Column 3: BP score
4.	CC semantic similarity score file (tab-delimited; no header):
  - Column 1: Protein 1
  - Column 2: Protein 2
  - Column 3: CC score
5.	Reference genome file (tab-delimited; no header):
  - Column 1: Entry (Accession)
  - Column 2: Entry name (ID)
  - Column 3: Gene names
  - Column 4: Protein names
  - Column 5: Organism ID
  - Column 6: Status
6.	Modules file (tab-delimited; with header):
  - Column 1: Module ID
  - Column 2: Number of module protein
  - Column 3: Module protein (Uniprot AC)

# Basic Usage

Step 1: Protein-protein interaction network (PPIN) construction

The example dataset is stored in the Dataset folder, all of the example outputs are provided in the PPIN folder.

Step 1-1: Calculate the R2 of PPIN under different PCC thresholds.

```bash
python3 1-1_R2_PPIN.py –e ./Dataset/9606_expPPI_network_OnlyAC.txt –p ./Dataset/COAD_PCC.txt –o ./PPIN/R2_COAD_PPIN.txt
```
- `h`: Get help with the commands.
- `e`: Experimental PPI file. 
- `p`: Gene co-expression file. 
- `o`: Output file for the calculated R2 values across all PCC thresholds.
