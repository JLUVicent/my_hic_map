#!/usr/bin/env python
# coding: utf-8
import numpy as np
import scipy.sparse as scisp
from math import log,exp,sqrt
import logging
import os
from raw_contact import ContactMatrix, ContactMatrix_LC

# package logger
logger = logging.getLogger(__name__)


class HiCzinMap:
    def __init__(self, path , contig_info , seq_map , norm_result , min_signal):
        '''
        perc: threshold of spurious contacts
        min_signal: minimum signal of acceptable contigs
        '''
        self.path = path
        self.seq_map = seq_map
        self.norm_result = norm_result
        self.signal = min_signal
        self.name = []
        self.site = []
        self.len = []
        self.cov = []
        self.tax = []

        for i in range(len(contig_info)):
            temp = contig_info[i]
            self.name.append(temp.name)
            self.site.append(temp.sites)
            self.len.append(temp.length)
            self.cov.append(temp.cov)
            self.tax.append(temp.tax)

        del contig_info

        self.name = np.array(self.name)
        self.site = np.array(self.site)
        self.len = np.array(self.len)
        self.cov = np.array(self.cov)
        self.tax= np.array(self.tax)
        
        self.cov[self.cov==0] = np.min(self.cov[self.cov!=0])
        
        self.norm()
        del self.site, self.cov

    def norm(self):
        self.seq_map = self.seq_map.tocoo()
        _map_row = self.seq_map.row
        _map_col = self.seq_map.col
        _map_data = self.seq_map.data
        _map_coor = list(zip(_map_row , _map_col , _map_data))
        coeff = self.norm_result[0:4]

        self.seq_map = self.seq_map.tolil()
        self.seq_map = self.seq_map.astype(np.float)
        for x,y,d in _map_coor:

            s1 = self.site[x]
            if s1 == 0:
                s1 = 1
            s2 = self.site[y]
            if s2 == 0:
                s2 = 1
            s = (log(s1*s2)-self.norm_result[5])/self.norm_result[6]

            l1 =  self.len[x]
            l2 =  self.len[y]
            l = (log(l1*l2)-self.norm_result[7])/self.norm_result[8]

            c1 = self.cov[x]
            c2 = self.cov[y]
            c = (log(c1*c2)-self.norm_result[9])/self.norm_result[10]
            
            d_norm = d/exp(coeff[0] + coeff[1]  * s  + coeff[2] * l + coeff[3]* c)
    
            if d_norm > self.norm_result[4]:
                self.seq_map[x , y] = d_norm
            else:
                self.seq_map[x , y] = 0
        del _map_row, _map_col, _map_data, _map_coor





class HiCzinMap_LC:
    def __init__(self, path , contig_info , seq_map , norm_result , min_signal):
        '''
        perc: threshold of spurious contacts
        min_signal: minimum signal of acceptable contigs 可接受的contigs的最小信号
        '''
        self.path = path
        self.seq_map = seq_map
        self.norm_result = norm_result
        self.signal = min_signal
        self.name = []
        self.len = []
        self.cov = []
        self.tax = []

        for i in range(len(contig_info)):
            temp = contig_info[i]#contig_info就是seq_info
            self.name.append(temp.name)
            self.len.append(temp.length)
            self.cov.append(temp.cov)
            self.tax.append(temp.tax)

        del contig_info

        self.name = np.array(self.name)#转换为数组
        self.len = np.array(self.len)
        self.cov = np.array(self.cov)
        self.tax= np.array(self.tax)

        
        self.norm()
        del self.cov



    def norm(self):
        self.seq_map = self.seq_map.tocoo()#稀疏矩阵
        _map_row = self.seq_map.row
        _map_col = self.seq_map.col
        _map_data = self.seq_map.data
        _map_coor = list(zip(_map_row , _map_col , _map_data))
        coeff = self.norm_result[0:4]

        self.seq_map = self.seq_map.tolil()#转换矩阵为list格式
        self.seq_map = self.seq_map.astype(np.float)#转换数据类型为浮点型
        for x,y,d in _map_coor:

            l1 =  self.len[x]
            l2 =  self.len[y]
            l = (log(l1*l2)-self.norm_result[4])/self.norm_result[5]

            c1 = self.cov[x]
            c2 = self.cov[y]
            c = (log(c1*c2)-self.norm_result[6])/self.norm_result[7]
            
            d_norm = d/exp(coeff[0] + coeff[1] * l + coeff[2]* c)
    
            if d_norm > self.norm_result[3]:
                self.seq_map[x , y] = d_norm
            else:
                self.seq_map[x , y] = 0
        del _map_row, _map_col, _map_data, _map_coor

if __name__=='__main__':
    from rpy2 import robjects
    r=robjects.r
    r.source('HiCzin.R')
    contig_file="/media/ubuntu/conda/vicent/HiCzin/my_614/contig_info.csv"
    valid_file="/media/ubuntu/conda/vicent/HiCzin/my_614/valid_contact.csv"
    #使用默认值
    thres=0.05
    #获得归一化结果
    norm_result=r.HiCzin(contig_file,valid_file,thres)
    print(type(norm_result))
    print("输出完成")
    #获取连接图脚本
    FASTA = '/media/ubuntu/conda/vicent/HiCzin/my_data_617/final.contigs.fa'
    BAM = '/media/ubuntu/conda/vicent/HiCzin/my_data_617/MAP_SORTED.bam'
    TAX = '/media/ubuntu/conda/vicent/HiCzin/my_data_617/vicent616_blast_out.csv'
    COV = '/media/ubuntu/conda/vicent/HiCzin/my_data_617/coverage616.txt'
    logger.info('Begin to test the contact map construction section...')
    runtime_defaults = {
        'min_len': 1000,
        'min_signal': 2,
        'min_mapq': 30,
        'min_match': 30,
        'thres': 0.05
    }
    # bam_file,  seq_file, tax_file, coverage_file , path , min_mapq=30, min_len=1000, min_match=30, min_signal=2
    cm = ContactMatrix_LC(BAM,
                    FASTA,
                    TAX,
                    COV,
                    './my_614',
                    min_mapq=runtime_defaults['min_mapq'],
                    min_len=runtime_defaults['min_len'],
                    min_match=runtime_defaults['min_match'],
                    min_signal=0)
    
    #归一化结果
    from utils import load_object, save_object, make_dir
    hzmap = HiCzinMap_LC('./my_614',
                cm.seq_info,
                cm.seq_map,
                norm_result,
                runtime_defaults['min_signal'])
    scisp.save_npz(os.path.join('/media/ubuntu/abc/vicent', 'Normalized_contact_matrix617.npz'), hzmap.seq_map.tocsr())
    save_object(os.path.join('/media/ubuntu/abc/vicent', 'HiCzin_normalized_contact617'), hzmap)
