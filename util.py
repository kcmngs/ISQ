#-*- conding:utf-8 -*-
# cython: language_level=3

import os
import glob
import re
import time
from time import strftime,localtime
import configparser


DEBUG = 0  # 在需要分析时效性的时候将该量置为1，否则置为0

#########################################################################
# 日志级别
#########################################################################
def log_info(level, info):
    """
    输入格式：
    level 日志级别，可以为：INFO WARN ERROR
    info 显示的内容
    返回：
    level YYYY-MM-dd H-M-S info
    """
    return "{level}\t{times}\t{info}".format(
        level=level,
        times=strftime("%Y-%m-%d %H:%M:%S", localtime()),
        info=info,
        )

