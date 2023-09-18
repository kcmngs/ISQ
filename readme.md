#### Parameter description
usage: IS_count.py

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT              ic result
  --anno ANNO           anno prefix
  --rmdup               rmdup[off]
  -c CONFIG, --config CONFIG
                        the config file
  -o OUTPUT             output prefix
  --log LOG             log

- -i Enter ic internal standard comparison results, reference format:test/294-BA-2-1.sample.ic
- --anno Enter the prefix of the file to be quantized，eg:test/294-BA-2-1.The program automatically finds the corresponding microbiology file *.anno.tsv
- --rmdup The default is not to go duplication,with this parameter use goduplication
- --config configuration file.The main configuration ic_group is used for each internal reference group setting，ic_group_nd is used to set the concentration gradient of the internal reference
- -o  Setting the output file prefix
- --log log file

#### Running Instructions
consultation:test/test.sh
