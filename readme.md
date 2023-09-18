#### 参数说明
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

- -i 输入ic内标比对结果，参考格式test/294-BA-2-1.sample.ic
- --anno 输入要定量的文件前缀，如test/294-BA-2-1，程序自动查找对应的微生物文件*.anno.tsv
- --rmdup 默认不进行去duplication,带上该参数时使用去duplication
- --config 配置文件，主要配置ic_group用于各个内参组设置，ic_group_nd用于设置内参浓度梯度
- -o  设置输出文件前缀
- --log 日志文件

#### 运行说明
参考test/test.sh
