# coding:utf-8
import os
from itertools import chain
import numpy as np
import xlsxwriter
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA,KernelPCA
from multiprocessing import Pool
import math
from collections import Counter
from scipy.optimize import linear_sum_assignment
from itertools import chain
import sys
from sklearn.metrics import accuracy_score
import json
sys.path.append("..")
from sklearn.cluster import KMeans
from sklearn import metrics
import gc

def g_con(matrix, chr_name):
    matrix1 = matrix
    for i in range(len(matrix)):
        cc = []
        index = []
        cc.append(matrix[i, :])
        for m in range(i - 1, i + 2):
            if m < 0 or m > len(matrix)-1 or m == i:
                continue
            cc.append(matrix[m, :])
            index.append(m)
        for n in range(len(matrix[0])):
            if n in index:
                continue
            if matrix[i, n] != 0 and i != n:
                cc.append(matrix[n, :])
        matrix1[i, :] = g_Mean(cc)
    m = np.triu(matrix1, 1)
    n = np.triu(matrix1, 1).T
    res = m + n + np.diag(np.diag(matrix1))
    return res

def g_Mean(List):
    mean = np.sum(np.array(List), axis=0)/len(List)
    return mean

def read_matrix(file_path):
    file = open(file_path)
    lines = file.readlines()
    a = []
    for line in lines:
        a.append(line.split())
    a = np.array(a).astype(float)
    return a

def matrix_list(matrix):
    return list(chain.from_iterable(matrix))


def original_select(matrix, prct):
    if prct > -1:
        thres = np.percentile(matrix, 100 - prct, axis=1)
        print(len(thres))#取前20%
        Q_concat = (matrix > thres[:, None])
    return Q_concat


def random_walk_imp(matrix, rp):
    row, _ = matrix.shape
    row_sum = np.sum(matrix, axis=1)
    for i in range(row_sum.shape[0]):
        if row_sum[i] == 0:
            row_sum[i] = 0.001
    nor_matrix = np.divide(matrix.T, row_sum).T
    Q = np.eye(row)
    I = np.eye(row)
    for i in range(30):
        Q_new = (1 - rp) * np.dot(Q, nor_matrix) + rp *I
        delta = np.linalg.norm(Q - Q_new)
        Q = Q_new.copy()
        if delta < 1e-6:
            break
    return Q

def processing_label(label, truth, n_clusters):
    point = 0
    k = 0
    label_order = []
    label_matrix = np.zeros((n_clusters, n_clusters), np.int)
    for i in range(len(truth)):
        if truth[point] == truth[i]:
            continue
        else:
            label_order.append(truth[point])
            counter = Counter(label[point:i])
            point = i
            for j in range(n_clusters):
                label_matrix[k][j] = counter[j]
            k += 1

    label_order.append(truth[point])
    counter = Counter(label[point:len(label) - 1])
    for j in range(n_clusters):
        label_matrix[n_clusters - 1][j] = counter[j]

    row_ind, col_ind = linear_sum_assignment(-label_matrix)
    for j in range(len(label)):
        for i in range(n_clusters):
            if label[j] == col_ind[i]:
                label[j] = label_order[i] + n_clusters
    return label - n_clusters


def con_ran(args):
    cell_id, type, chr_name, index, rp = args
    file_path = "./contact_626/%s/cell_%s_%s.txt" % (type, str(cell_id), chr_name)
    chr_file = open(file_path)
    lines = chr_file.readlines()
    contact_matrix = np.zeros((index, index))
    for line in lines:
        bin1, bin2, num = line.split()
        contact_matrix[int(bin1), int(bin2)] += int(num)
        if bin1 != bin2:
            contact_matrix[int(bin2), int(bin1)] += int(num)
    imp_matrix1 = g_con(contact_matrix, chr_name)
    r_matrix1 = random_walk_imp(imp_matrix1, rp)
    r_path1 = "./temp/Adjacent/Adj_rrw_matrix/%s/cell_%s_%s.txt" % (type, str(cell_id), chr_name)
    np.savetxt(r_path1, r_matrix1, fmt='%f', delimiter=' ')

def del_file(path):
    ls = os.listdir(path)
    for i in ls:
        c_path = os.path.join(path, i)
        if os.path.isdir(c_path):
            del_file(c_path)
        else:
            os.remove(c_path)

def dele_temp():
    del_file('./temp/Adjacent/Adj_rrw_matrix/GM12878')
    del_file('./temp/Adjacent/Adj_rrw_matrix/K562')
    del_file('./temp/Adjacent/Adj_rrw_matrix/HeLa')
    del_file('./temp/Adjacent/Adj_rrw_matrix/HAP1')


def main():
    cell_num_path = open('./cell_num.json', 'r')
    types = json.load(cell_num_path)
    rp = 0.9
    f = open("./combo_hg19.genomesize")
    index = {}
    lines = f.readlines()
    for line in lines:
        chr_name, length = line.split()
        chr_name = chr_name.split('_')[1]
        max_len = int(int(length)/1000000)
        index[chr_name] = max_len + 1
    f.seek(0, 0)
    print(index)
    f.close()
    #########################################################################################################领近平滑+重启随机游走
    p = Pool(10)
    for type in types:
        for c_num in range(types[type]):
            cell_id = c_num + 1
            args = [[cell_id, type, chr_name, index[chr_name], rp] for chr_name in index]# 进程池，最大的进程池为10
            print(args)
            p.map(con_ran, args)# 调用impute_cpu函数，result包含[细胞，矩阵]
    p.close()
    gc.collect()
    # #########################################################################################################一次KPCA 用sigmoid和linear -gc
    kernel = "linear"
    cell_matrix = []
    new_chr_matrix = []
    for chr_name in index:
        new_chr_matrix.clear()
        pca_index = []
        for type in types:
            for c_num in range(types[type]):
                cell_id = c_num + 1
                r_path = "./temp/Adjacent/Adj_rrw_matrix/%s/cell_%s_%s.txt" % (type, str(cell_id), chr_name)
                r_matrix = read_matrix(r_path)
                new_chr_matrix.append(matrix_list(r_matrix))
                cindex = (type, cell_id)
                pca_index.append(cindex)
        new_chr_matrix = np.array(new_chr_matrix)
        cell_matrix.append(original_select(new_chr_matrix, 13))
        new_chr_matrix = new_chr_matrix.tolist()
    dele_temp()
    cell_matrix = np.concatenate(cell_matrix, axis=1)
    pca = KernelPCA(n_components=min(cell_matrix.shape) - 1, kernel=kernel) #min(C_matrix.shape) - 1 # laplacian "linear" | "poly" | "rbf" | "sigmoid" | "cosine" | "precomputed"
    pca.fit(cell_matrix)
    pca_matrix = pca.fit_transform(cell_matrix)
    label = []
    for i in pca_index:
        if i[0] == 'GM12878':
           label.append(0)
        elif i[0] == 'HAP1':
            label.append(1)
        elif i[0] == 'HeLa':
            label.append(2)
        elif i[0] == 'K562':
            label.append(3)


if __name__ == "__main__":
    main()
