#-*- conding:utf-8 -*-
# cython: language_level=3

import sys
import codecs
import numpy as np
import math
import statsmodels.api as sm
import pandas as pd
#import matplotlib.pyplot as plt
import matplotlib
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import scipy.stats as sst
import configparser
import argparse
import os
from util import log_info
matplotlib.use('agg')

class ConfigParserUper(configparser.ConfigParser):
    def optionxform(self,optionstr):
        return optionstr

def ols_formula(X, Y, picture=None, log=None):
    n = len(X)
    x_mean = np.mean(X)
    lxx = np.sum((X-x_mean)**2)
    y_mean = np.mean(Y)
    lyy = np.sum((Y-y_mean)**2)
    lxy = np.sum((X-x_mean)*(Y-y_mean))

    # X列加入常数项1
    X_n = sm.add_constant(X)
    # Y=因变量 X_n是自变量
    model = sm.OLS(Y,X_n)
    # 拟合结果
    results = model.fit()
    # 打印拟合摘要
    if log:
        with open(log, 'w') as f:
            f.write(results.summary().as_text())
    # 回归系数
    params = results.params
    a = params[0]
    b = params[1]
    p_values = results.pvalues
    Qe = lyy - b * lxy
    # 总体方差
    sigma_est2 = Qe / (n - 2)
    sigma_est = np.sqrt(sigma_est2)
    test_level = sst.t.ppf(1 - 0.05/2, df=n-2)
    if picture:
        # 标准差 置性区间
        prstd, ivLow, ivUp = wls_prediction_std(results)
        fig, ax = matplotlib.pyplot.subplots(figsize=(10, 8))
        ax.plot(X, Y, 'o', label="data")
        ax.plot(X, results.fittedvalues, 'r-', label="OLS")
        ax.plot(X, ivUp, '--',color='blue',label="95%upConf")
        ax.plot(X, ivLow, '--',color='orange',label="95%lowConf")
        matplotlib.pyplot.xlabel('ris')
        matplotlib.pyplot.ylabel('cis')
        R2="R2="+str('%.4f' % results.rsquared)
        p_a="P(a)_values="+str('%.4f' % p_values[0])
        p_b="P(b)_values="+str('%.4f' % p_values[1])
        y="y="+str('%.4f' % a)+"+"+str('%.4f' % b)+"x"
        matplotlib.pyplot.text(2,-1, R2)
        matplotlib.pyplot.text(2,-1.2, p_a)
        matplotlib.pyplot.text(2,-1.4, p_b)
        matplotlib.pyplot.text(2,-1.6,y)
        ax.legend(loc='lower right')
        matplotlib.pyplot.title('OLS linear regression ')
        matplotlib.pyplot.savefig(picture)
    return ({
        "a":a,
        "b":b,
        "p_values": p_values,
        "sigma_est": sigma_est,
        "test_level": test_level,
        "lxx": lxx,
        "n": n,
        "x_mean": x_mean,
        "R2": results.rsquared,
        })


def IC_ols_formula(ic, cf, output, rmdup=True, log=None):
    '''
    计算线性共识
    输入:
        内标序列比对结果
            格式：
                sample          ic_num  map_reads_rmdup map_reads std_map_reads_rmdup std_map_reads
                294-BA-2-1      IKC_f1  250922  250922  199081  199081
        config
            [ic_group]
            f1_f3_f4=IKC_f1,IKC_f3,IKC_f4
            f6_f8_f9=IKC_f6,IKC_f8,IKC_f9
            f10_f11_f14=IKC_f10,IKC_f11,IKC_f14
            f15_f18_f23=IKC_f15,IKC_f18,IKC_f23
            f24_f27_f28=IKC_f24,IKC_f27,IKC_f28
            f1_f3_f4_nd=0.05
            f6_f8_f9_nd=0.01
            f10_f11_f14_nd=0.002
            f15_f18_f23_nd=0.0004
            f24_f27_f28_nd=0.00008
    输出:None 或 dict
    '''
    cf = ConfigParserUper()
    cf.read(args.config, encoding='utf-8')
    IC_group_list = []
    IC_group_name = []
    IC_group_list_num = []
    for ic_k,ic_v in cf.items('ic_group'):
        IC_group_name.append(ic_k)
        IC_group_list.append(ic_v.split(','))
        IC_group_list_num.append([0,0])
    with codecs.open(ic, 'r', 'gbk') as f:
        for line in f:
            line = line.strip().split("\t")
            for i in range(len(IC_group_list)):
                if line[1] in IC_group_list[i]:
                    #map reads
                    IC_group_list_num[i][0] += int(line[3])
                    #rmdup map reads
                    IC_group_list_num[i][1] += int(line[2])

    # to log(10)
    TR = 2
    IC_group_log = []
    IC_group_nd_log = []
    total_ic_reads = 0
    for i in range(len(IC_group_name)):
        ic_name = IC_group_name[i]
        if rmdup:
            index = 1
        else:
            index = 0
        IC_group_nd_log.append(math.log(cf.getfloat('ic_group_nd', ic_name+"_nd") * TR , 10))
        IC_group_log.append( 0 if IC_group_list_num[i][index] <= 0 else math.log(IC_group_list_num[i][index], 10))
        total_ic_reads += IC_group_log[-1]

    # 无内标结果
    if total_ic_reads <= 0:
        return None
    return ols_formula(X=np.array(IC_group_log), Y=np.array(IC_group_nd_log), picture=output+".ols.png", log=log)


def static_copies(genome_len, reads, ic_parm):
    if ic_parm is None or 'a' not in ic_parm or 'b' not in ic_parm:
        return('-', '-','-')
    if reads <= 0:
        return(0,0,0)
    #zl = genome_len * (10**-6)
    # 预测浓度
    Sy_log10 = ic_parm['a'] + ic_parm['b'] * math.log(reads ,10)
    # 种求指数：理论投入量
    Sylgtouru=pow(10,Sy_log10) * pow(10, 9)
    # 投入量换算拷贝数 并且向上取整
    beichushu = genome_len * 1.022 * pow(10, -3) 
    #Skb = math.ceil((Sylgtouru/(zl*(10**6)*(1.022*10**-9)*(10**-3)))
    Skb = Sylgtouru / ( genome_len * 1.022 * pow(10, -3))
    Skb = math.ceil(Skb)
    # 95置性区间
    Sytouru_down = Sy_log10 - ic_parm['test_level'] * ic_parm['sigma_est'] * np.sqrt(1 + 1 / ic_parm['n'] + ((math.log(reads,10) - ic_parm['x_mean']) ** 2 / ic_parm['lxx']))
    Sylgtouru_down = pow(10,Sytouru_down) * pow(10, 9)
    #Spi_down = Sylgtouru_down / (zl * (10**6) * (1.022 * 10**-9) * (10**-3))
    Spi_down = Sylgtouru_down / (genome_len * 1.022 * pow(10, -3))
    Spi_down = math.ceil(Spi_down)
    Sytouru_up = Sy_log10 + ic_parm['test_level'] * ic_parm['sigma_est'] * np.sqrt(1 + 1 / ic_parm['n'] + (( math.log(reads,10) - ic_parm['x_mean']) ** 2 / ic_parm['lxx']))
    Sylgtouru_up = pow(10,Sytouru_up) * pow(10, 9)
    #Spi_up = Sylgtouru_up / (zl * (10**6) * (1.022 * 10**-9) * (10**-3))
    Spi_up = Sylgtouru_up / (genome_len * 1.022 * pow(10, -3))
    Spi_up = math.ceil(Spi_up)
    return(Skb, Spi_up,Spi_down)

def target_gene(ic_parm,input_file,output_file):
    title = None
    linenum = 0
    with codecs.open(input_file, mode='r', encoding='gbk') as fin, open(output_file, 'w', encoding='gbk') as fout:
        for line in fin:
            linenum += 1
            line = line.strip().split("\t")
            if linenum == 1:
                title = line
                title.extend(["species_copy", "species_conf","genus_copy", "genus_conf"])
                fout.write("\t".join(title)+"\n")
                continue
            # 建立位置索引
            species_idx = title.index("species_reads")
            genus_idx = title.index("genus_reads")
            coverage_idx = title.index("coverage")
            # 种reads数据
            species_reads = 0 if line[species_idx] == '-' else int(line[species_idx])
            # 属reads数据
            gennus_reads = 0 if line[genus_idx] == '-' else int(line[genus_idx])
            # 基因长度数据
            genome_len = int(line[coverage_idx].strip().split("/")[1])
            # 种拷贝数， 95置信区间
            Skb, Spi_up, Spi_down = static_copies(genome_len=genome_len, reads=species_reads, ic_parm=ic_parm)
            Gkb, Gpi_up, Gpi_down = static_copies(genome_len=genome_len, reads=gennus_reads, ic_parm=ic_parm)
            fout.write("{lines}\t{Skb}\t{Spi_down};{Spi_up}\t{Gkb}\t{Gpi_down};{Gpi_up}\n".format(
                lines="\t".join(line),
                Skb=Skb,
                Spi_down=Spi_down,
                Spi_up=Spi_up,
                Gkb=Gkb,
                Gpi_down=Gpi_down,
                Gpi_up=Gpi_up,
                )
            )
                

def virus_gene(ic_parm, input_file,output_file):
    title = None
    linenum = 0
    with codecs.open(input_file,mode='r',encoding='gbk') as fin, open(output_file, 'w', encoding='gbk') as fout:
        for line in fin:
            linenum += 1
            line = line.strip().split("\t")
            if linenum == 1:
                title = line
                title.extend(["subtype_copy", "subtype_conf", "serotype_copy", "serotype_conf", "species_copy", "species_conf"])
                fout.write("\t".join(title)+"\n")
                continue
            # 建立位置索引
            subtype_idx = title.index("subtype_reads")
            serotype_idx = title.index("serotype_reads")
            species_idx = title.index("species_reads")
            coverage_idx = title.index("coverage")
            # 亚型序列数
            subtype_reads = 0 if line[subtype_idx] == '-' else int(line[subtype_idx])
            # 型序列数
            serotype_reads = 0 if line[serotype_idx] == '-' else int(line[serotype_idx])
            # 种序列数
            species_reads =  0 if line[species_idx] == '-' else int(line[species_idx])
            # 基因组长度
            genome_len = int(line[coverage_idx].strip().split("/")[1])
            subtype_kb, subtype_pi_up, subtype_pi_down = static_copies(genome_len=genome_len, reads=subtype_reads, ic_parm=ic_parm)
            serotype_kb, serotype_pi_up, serotype_pi_down = static_copies(genome_len=genome_len, reads=serotype_reads, ic_parm=ic_parm)
            Skb, Spi_up, Spi_down = static_copies(genome_len=genome_len, reads=species_reads, ic_parm=ic_parm)
            fout.write("{lines}\t{subtype_kb}\t{subtype_pi_down};{subtype_pi_up}\t{serotype_kb}\t{serotype_pi_down};{serotype_pi_up}\t{Skb}\t{Spi_down};{Spi_up}\n".format(
                lines="\t".join(line),
                subtype_kb=subtype_kb,
                subtype_pi_down=subtype_pi_down, 
                subtype_pi_up=subtype_pi_up,
                serotype_kb=serotype_kb,
                serotype_pi_down=serotype_pi_down,
                serotype_pi_up=serotype_pi_up,
                Skb=Skb,
                Spi_down=Spi_down,
                Spi_up=Spi_up,
                )
            )
if __name__=='__main__':
    parser = argparse.ArgumentParser(usage="%(prog)s")
    parser.add_argument('-i', required=True, help='ic result', dest='input')
    parser.add_argument('--anno', required=True, help='anno prefix', dest='anno')
    parser.add_argument('--rmdup', required=False, help='rmdup[off]', dest='rmdup', action='store_true', default=False)
    parser.add_argument('-c', '--config', action='store', required=True, help='the config file')
    parser.add_argument('-o', required=True, help='output prefix', dest='output')
    parser.add_argument('--log',required=False, help='log', dest='log')
    args = parser.parse_args()
    if not os.path.exists(args.input):
        sys.stderr.write(log_info(level="WARN", info="%s not exists" % args.input)+"\n")
        sys.exit(1)
    if not os.path.exists(args.config):
        sys.stderr.write(log_info(level="WARN", info="%s not exists\n" % args.config)+"\n")
        sys.exit(1)
    result = IC_ols_formula(ic=args.input, cf=args.config, rmdup=args.rmdup, output=args.output, log=args.log)
    with open(args.output+".result", 'w') as f:
        if result is not None:
            f.write("a\t{a}\nb\t{b}\nR2\t{R2}\na_pvalue\t{a_pvalue}\nb_pvalue\t{b_pvalue}\nQC\t{result}\n".format(
                a="%.4f" % result['a'], 
                b="%.4f" % result['b'], 
                R2="%.4f" % result['R2'], 
                a_pvalue="%.4f" % result['p_values'][0], 
                b_pvalue="%.4f" % result['p_values'][1],
                result="Y" if result['R2'] >= 0.95 and result['p_values'][0] <= 0.05 and result['p_values'][1] <= 0.05 else "N",
                )
            )
        else:
            f.write("a\t-\nb\t-\nR2\t-\na_pvalue\t-\nb_pvalue\t-\nQC\tN\n")
                
    # anno
    if os.path.exists(args.anno + ".Bacteria.anno.tsv"):
        target_gene(
            ic_parm=result,
            input_file=args.anno+".Bacteria.anno.tsv",
            output_file=args.output+".Bacteria.anno.IS.tsv"
            )
    if os.path.exists(args.anno + ".Fungi.anno.tsv"):
        target_gene(
            ic_parm=result,
            input_file=args.anno+".Fungi.anno.tsv",
            output_file=args.output+".Fungi.anno.IS.tsv",
            )
    if os.path.exists(args.anno + ".Parasite.anno.tsv"):
        target_gene(
            ic_parm=result,
            input_file=args.anno+".Parasite.anno.tsv",
            output_file=args.output+".Parasite.anno.IS.tsv"
            )
    if os.path.exists(args.anno + ".Virus.anno.tsv"):
        virus_gene(
            ic_parm=result,
            input_file=args.anno+".Virus.anno.tsv",
            output_file=args.output+".Virus.anno.IS.tsv"
            )
    # filter raw
    if os.path.exists(args.anno + ".Bacteria.anno.raw.filter.tsv"):
        target_gene(
            ic_parm=result,
            input_file=args.anno+".Bacteria.anno.raw.filter.tsv",
            output_file=args.output+".Bacteria.anno.raw.filter.IS.tsv"
            )
    if os.path.exists(args.anno + ".Fungi.anno.raw.filter.tsv"):
        target_gene(
            ic_parm=result,
            input_file=args.anno+".Fungi.anno.raw.filter.tsv",
            output_file=args.anno+".Fungi.anno.raw.filter.IS.tsv"
            )
    if os.path.exists(args.anno + ".Parasite.anno.raw.filter.tsv"):
        target_gene(
            ic_parm=result,
            input_file=args.anno+".Parasite.anno.raw.filter.tsv",
            output_file=args.anno+".Parasite.anno.raw.filter.IS.tsv"
            )
    if os.path.exists(args.anno + ".Virus.anno.raw.filter.tsv"):
        virus_gene(
            ic_parm=result,
            input_file=args.anno+".Virus.anno.raw.filter.tsv",
            output_file=args.anno+".Virus.anno.raw.filter.IS.tsv",
            )
    # raw
    if os.path.exists(args.anno + ".Bacteria.anno.raw.tsv"):
        target_gene(
            ic_parm=result,
            input_file=args.anno+".Bacteria.anno.raw.tsv",
            output_file=args.output+".Bacteria.anno.raw.IS.tsv"
            )
    if os.path.exists(args.anno + ".Fungi.anno.raw.tsv"):
        target_gene(
            ic_parm=result,
            input_file=args.anno+".Fungi.anno.raw.tsv",
            output_file=args.anno+".Fungi.anno.raw.IS.tsv"
            )
    if os.path.exists(args.anno + ".Parasite.anno.raw.tsv"):
        target_gene(
            ic_parm=result,
            input_file=args.anno+".Parasite.anno.raw.tsv",
            output_file=args.anno+".Parasite.anno.raw.IS.tsv"
            )
    if os.path.exists(args.anno + ".Virus.anno.raw.tsv"):
        virus_gene(
            ic_parm=result,
            input_file=args.anno+".Virus.anno.raw.tsv",
            output_file=args.anno+".Virus.anno.raw.IS.tsv",
            )
