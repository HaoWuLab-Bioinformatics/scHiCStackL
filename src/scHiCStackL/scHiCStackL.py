# coding:utf-8
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from itertools import chain
import numpy as np
from sklearn.model_selection import  KFold
from sklearn.linear_model import RidgeClassifier
from sklearn.naive_bayes import GaussianNB
from scipy import stats
import json
import time

def read_matrix(file_path):
    file = open(file_path)
    lines = file.readlines()
    a = []
    for line in lines:
        a.append(line.split())
    a = np.array(a).astype(float)
    return a

def vote(temp_test, valid_length):
    New_Test_Data = np.zeros((valid_length, 2))
    for j in range(valid_length):
        for k in range(2):
            pr = []
            for l in temp_test.keys():
                pr.append(temp_test[l][j, k])
            New_Test_Data[j, k] = stats.mode(pr)[0][0]
    return New_Test_Data

def original_select(matrix, prct):
    if prct > -1:
        thres = np.percentile(matrix, 100 - prct, axis=1)
        print(len(thres))
        Q_concat = (matrix > thres[:, None])
    return Q_concat

def matrix_list(matrix):
    return list(chain.from_iterable(matrix))


def LR1(X_train_train, y_train_train,X_train_test,valid_X, C):
    clf = LogisticRegression(max_iter=2500, C=C)
    clf.fit(X_train_train, y_train_train)
    train_label_pred = clf.predict(X_train_test)
    test_label_pred = clf.predict(valid_X)
    return train_label_pred, test_label_pred

def Ridge(X_train_train, y_train_train,X_train_test,valid_X, al):
    clf = RidgeClassifier(alpha=al, max_iter=2500)
    clf.fit(X_train_train, y_train_train)
    train_label_pred = clf.predict(X_train_test)
    test_label_pred = clf.predict(valid_X)

    return train_label_pred, test_label_pred

def read_json(path, ndim):
    f = open(path, 'r')
    parameters = json.load(f)[str(ndim)]
    return np.array(parameters)


def predict(Data,Label, cell_num, ndim):
    """
    Two-layer stacking ensemble model
    train and predict data
    :param Data: array-like, shape=(n_samples, n_features)
    :param Label: array-like, shape=(n_samples, 1)
    :param cell_num:  int, 626 (cell_type = “human”) | 2655 (cell_type = “human”) | 178(cell_type = “mouse”)
    :param ndim: int, [10,49]
                ndim.Default=40.
    :return:  ARI, Acc, MCC, F1, Precision, NMI
    """
    X = np.array(Data)
    Y = (np.array(Label)).astype(np.uint8)
    n_splits1 = 5
    n_splits2 = 10
    folds = KFold(n_splits1, shuffle=True, random_state=0).split(X, Y)
    C_matrix = read_json("./parameters/"+str(cell_num)+"/C.json", ndim)
    al_matrix = read_json("./parameters/" + str(cell_num) + "/al.json", ndim)
    ARI = Acc = NMI = F1 = Precision = MCC = 0
    temp_test = {}
    for predict_num, (trained, valided) in enumerate(folds):
        train_y, train_X = Y[trained], X[trained]
        valid_y, valid_X = Y[valided], X[valided]
        New_Train_Data = []
        Train_Label = []
        kf = KFold(n_splits=n_splits2, shuffle=True, random_state=0)
        for iter_num, (train_index, test_index) in enumerate(kf.split(train_X, train_y)):
            X_train_train, X_train_test = train_X[train_index], train_X[test_index]
            y_train_train, y_train_test = train_y[train_index], train_y[test_index]
            C = C_matrix[predict_num, iter_num]
            al = al_matrix[predict_num, iter_num]
            train_label_pred1, test_label_pred1 = LR1(X_train_train, y_train_train, X_train_test, valid_X, C)
            train_label_pred2, test_label_pred2 = Ridge(X_train_train, y_train_train, X_train_test, valid_X, al)
            temp_train = np.c_[train_label_pred1, train_label_pred2]
            temp_test[iter_num] = np.c_[test_label_pred1, test_label_pred2]
            New_Train_Data.append(temp_train)
            Train_Label.append(y_train_test)
        New_Train_Data = np.concatenate(New_Train_Data).astype(np.uint8)
        Train_Label = np.concatenate(Train_Label).astype(np.uint8)
        New_Test_Data = vote(temp_test, len(valid_y))
        Test_Label = np.array(valid_y).astype(np.uint8)
        ##################################################################################
        model = GaussianNB()
        model.fit(New_Train_Data, Train_Label)
        label_pred = model.predict(New_Test_Data)
        ##################################################################################
        Precision += metrics.precision_score(Test_Label, label_pred, average='weighted')
        MCC += metrics.matthews_corrcoef(Test_Label, label_pred)
        Acc += metrics.accuracy_score(Test_Label, label_pred)
        F1 += metrics.f1_score(Test_Label, label_pred, average='weighted')
        ARI += metrics.adjusted_rand_score(Test_Label, label_pred)
        NMI += metrics.normalized_mutual_info_score(Test_Label, label_pred)
    q = np.array([ARI,  Acc, MCC, F1, Precision, NMI])/5
    return q


def stackling_model(cell_type="human", cell_num = 626, ndim = 40):
    """
    Load PCA file and cell labels
    :param cell_type: “human” | “mouse”
            cell_type. Default="human".
    :param cell_num:  int, 626 (cell_type = “human”) | 2655 (cell_type = “human”) | 178(cell_type = “mouse”)
            cell_num. Default=626.
    :param ndim: int, [11,50]
                ndim.Default=40.
    :return: ARI, Acc, MCC, F1, Precision, NMI
    """
    Data_path = "./Files/PCA_file/%s/%s/pca_cell_file.txt" %(cell_type, str(cell_num))
    Data = read_matrix(Data_path)[:, :ndim]
    Lable_path = "./Files/PCA_file/%s/%s/label.txt" %(cell_type, str(cell_num))
    Label = read_matrix(Lable_path).flatten()
    Results = predict(Data,Label,cell_num,ndim-1)
    ARI  = Results[0]
    Acc = Results[1]
    MCC = Results[2]
    F1 = Results[3]
    Precision = Results[4]
    NMI = Results[5]
    return ARI, Acc, MCC, F1, Precision, NMI