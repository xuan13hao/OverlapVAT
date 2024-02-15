# Introduction
## Install
```console
mkdir build
cd build
cmake ..
make
```
## DNA
## DNA for multiple threads
```console
./VAT makevatdb --dbtype nucl --in ../data/strep_2x.fasta -d strep_2x 
/usr/bin/time -v ./VAT dna -d strep_2x.vatf -q ../data/strep_2x.fasta -p 1 > log
```

```
## Protein
```console
./VAT makevatdb --dbtype prot --in ../data/strep_2x.fasta -d strep_2x
/usr/bin/time -v ./VAT protein -d strep_2x.vatf -q ../data/strep_2x.fasta -p 1 > log
```
vim match
```

## DNA for minimap2
```console
/usr/bin/time -v ./VAT dna -r ../data/strep_2x.fasta -q ../data/strep_2x.fasta --mini -p 1 > log
```

