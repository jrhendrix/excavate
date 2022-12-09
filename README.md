# Excavate

A tool for filtering Megalodon methylation output into a GFF (General Feature Format) file. 

## Getting started
### Requirements
* Python 3.7+
* numpy
* matplotlib

### Installation
Download executable file from GitHub and run. No special installation required.

### Usage
As input, Excavate takes as input a Megalodon database file in tsv format.

Example1: basic
``` python excavate -i mega_db.tsv ```

Example 2: do not export gff file
``` python excavate -i mega_db.tsv -g ```

Example 3: Set probability threshold to 80%
``` python excavate -i mega_db.tsv -t 80 ```

### Output Files
* Modification stats
* Modification stats for significant only
* GFF of significant modifications
