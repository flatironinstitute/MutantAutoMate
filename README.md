# MutantAutoMate

- The script `getPDB.py` will find the isoform of a gene and mutant of interest and find and download the relevant PDB.
- The script `getPDF2.py` will generate a PDF with relevant descriptions of a given mutant.
- The script `snapshot.sh` will produce an image of your mutant structure.

Download dependencies after cloning repository with 

`pip install -e.`

Call `snapshot.sh` through the terminal by:

`./snapshot.sh 'file1' 'file2' `

To run a code in the src folder:

`python -m mutantautomate getPDB.py`

This code for NRXN1 for residue P89 will give the output(isoforms list):

P55196
Q62765
Q62889
Q8N2Q7
Q07666
Q9NZ94
Q99K10
Q91V33
Q60749
O35889
Q9QZQ1
P55196-5
P55196-1
P55196-2
P55196-3
P55196-6
Q62765-3
Q62765-2
Q62765-4
Q62889-2
Q62889-4
Q62889-3


## Examples

```bash
python ./final.py --gene-name NLGN1 --residue1 D --position 140 --residue2 Y --top-isoforms True
```

