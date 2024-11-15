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
docker compose run python bash
cd src/mutantautomate
python ./final.py --gene-name NLGN1 --residue1 D --position 140 --residue2 Y --top-isoforms True
```

## Frontend / Web Server

[app2.py](src/mutantautomate/app2.py) is a Flask app that serves a single index page, [index2.html](src/mutantautomate/templates/index.html).

The Flask app uses [`process.py`](src/mutantautomate/process.py) and [`pdb_helpers.py`](src/mutantautomate/pdb_helpers.py) to do the back-end processing. `process.py` has a `process` function which is a Python *generator*. The front-end opens up an `EventSource`, and the `process` function yields JSON events back to the browser.

The front-end uses [Preact](https://preactjs.com) and other JS libraries to render the front end. There is no build step, everything is rendered in-browser and fetched from the JS CDN called [esm.sh](https://esm.sh).