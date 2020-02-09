import numpy as np
from scipy import stats
import scipy.io
import os
import fun
import Tkinter, tkFileDialog
import random
import pywt
from sklearn import svm
from sklearn import tree
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.neural_network import MLPClassifier
from sklearn.externals import joblib

##################################################
## load recording files from selected directory ##
##################################################

# load answer sheet before the loop
filename_path = 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\data\\blvg_trace_exp2_all\\answer\\answer_sheet.mat'
mat_answer = scipy.io.loadmat(filename_path)
answer_sheet = mat_answer['answer_sheet']

# load all accuracy
all_accuracy = np.load('C:\\Users\\mz86\\Desktop\\fNIRS ML\\trained_classifiers\\accuracy.npy')
# load all loss
all_loss = np.load('C:\\Users\\mz86\\Desktop\\fNIRS ML\\trained_classifiers\\loss.npy')

# load data_matrix from each recording site
# for each file 
# divided to 2 parts (80% training, 20% testing)
# do this for n times
# each time 80% used to train SVM
# test this SVM on the 20% data
# output the accuracy each time and calculate average

file_path = os.path.dirname(os.path.realpath(__file__))
root = Tkinter.Tk()
root.withdraw()
dirname = tkFileDialog.askdirectory(parent=root,initialdir=file_path,title='Please select a directory')
filename_list=[]


four_sit_result = []
four_sit_result_DT = []
four_sit_result_nn_ITD= []
four_sit_result_nn_noITD= []



for filename in os.listdir(dirname):
    
    # load matlab file and save filename
    print(filename)
    filename_list.append(filename)
    site_name = filename.split('.')[0]
    
    filename_path=dirname+'/'+filename
    mat = scipy.io.loadmat(filename_path)
    data_matrix = mat[site_name]
    
    # load all trained classifier for this region and do the test
    path = 'C:\\\Users\\mz86\\Desktop\\fNIRS ML\\trained_classifiers\\' + site_name + '\\'
    
    # divide to 2 parts: training and testing
    # random divide for n times
    # each time to training and classification
    all_per_correct = []
    all_per_correct_DT = []
    all_per_correct_nn = []
    
    all_ones_nn = []
    
    # for each site chose top m (odd number) to form a voting system
    # highers accuracy or least loss
    m = 19
    if site_name == 'A_matrix':
        accuracy = all_accuracy[0]
        loss = all_loss[0]
    elif site_name == 'B_matrix':
        accuracy = all_accuracy[1]
        loss = all_loss[1]
    elif site_name == 'C_matrix':
        accuracy = all_accuracy[2]
        loss = all_loss[2]
    elif site_name == 'D_matrix':
        accuracy = all_accuracy[3]
        loss = all_loss[3]
    
    # find top m indexed
    selected_nn_index = accuracy.argsort()[-m:]
    selected_nn_index_loss = loss.argsort()[0:m]
    
    # only repeat m time with 100 as testing
    #repeat_num = 1
    for n in selected_nn_index:
        
        # how many rows, all of them is active listening
        row_num = np.shape(data_matrix)[0]  
        
        # file path of the trained classifier
        classifier_path = path + str(n)
           
        ## load 3 trained classifier
        clf_svm = joblib.load(classifier_path + '\\clf_svm.pkl')
        clf_DT = joblib.load(classifier_path + '\\clf_DT.pkl')
        clf_nn = joblib.load(classifier_path + '\\clf_nn.pkl')
        ## load selected feature index
        feature_index = np.load(classifier_path + '\\feature_index.npy')
        
        ## wavelet transform of each row of blvg matrix for this site
        ## select top 10 least likely to be random as feature later
        ## and build matrix for SVM classification
        wave_coef_matrix=[]
        for blvg_trace in data_matrix:            
            # wave_coef=fun.discreteHaarWaveletTransform(blvg_trace)       
            all_level_coef = pywt.wavedec(blvg_trace, 'haar', level=4, mode='per')  
            wave_coef = np.concatenate( all_level_coef, axis=0 )    
            wave_coef_matrix.append(wave_coef)        
        wave_coef_matrix=np.asarray(wave_coef_matrix)
        # delete last column that is all 0
        wave_coef_matrix = np.delete(wave_coef_matrix, -1, axis=1)
        
        ### to train on training set only
        #wave_coef_matrix_training = wave_coef_matrix[training_index,:]
        #wave_coef_matrix_testig = wave_coef_matrix[testing_index,:]
        
        ## from wave_coef_matrix choose 10 columns least likely to be random
        # using Kolmogorov-Smirnov (KS) test
        # http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.653.6292&rep=rep1&type=pdf
        all_feature_test_result = []
        for ii in range(0, wave_coef_matrix.shape[1]):
            feature = wave_coef_matrix[:,ii].astype(float)
            # run KS test on this feature
            # https://stackoverflow.com/questions/17901112/using-scipys-stats-kstest-module-for-goodness-of-fit-testing
            test_result = stats.kstest(feature, 'norm', args=(np.mean(feature), np.std(feature)))
            # test_result 0: test static 1: p value
            all_feature_test_result.append(test_result[0])
        all_feature_test_result = np.asarray(all_feature_test_result)
            
        ## pick the n largest p value column as selected features to build the matrix
        # selected_feature_index = all_feature_test_result.argsort()[-10:][::1]
        feature_num = 25
        selected_feature_index = all_feature_test_result.argsort()[:feature_num]
        selected_feature_index = feature_index;
    
        ## do PCA analysis after wavelet anaylis
        ## inital 8 PC
        ## plot top 2 PCA component in 2d graph 
        #pca_a = PCA(n_components=8)
        #pca_a.fit(wave_coef_matrix)
        #eigen_a=pca_a.explained_variance_ratio_
        #eigen_vector_a=pca_a.components_
        #blvg_trace_r_a = pca_a.fit(wave_coef_matrix).transform(wave_coef_matrix)
     
        # use clf to perdict 
        #svm
        out_put = clf_svm.predict(wave_coef_matrix[:,selected_feature_index])
        #DT
        out_put_DT = clf_DT.predict(wave_coef_matrix[:,selected_feature_index])
        #nn
        out_put_nn = clf_nn.predict(wave_coef_matrix[:,selected_feature_index])
        
        # save the output 0 and 1 
        all_ones_nn.append(out_put_nn)
        
    # for each column the majority of 1 or 0 wins the result  
    perdiction = []
    all_ones_nn = np.asarray(all_ones_nn)
    
    for jj in range(0, all_ones_nn.shape[1]):
        voting_results = np.count_nonzero(all_ones_nn[:,jj])
        if voting_results >= (m+1)/2:
            perdiction.append(1)
        else:
            perdiction.append(0)
                
    # calculate pertange correct
    perdiction = np.asarray(perdiction)
    
    per_correct_ITD = np.sum(perdiction[0:13].astype(float))/(row_num/2)
    per_correct_noITD = np.sum(perdiction[14:27].astype(float))/(row_num/2)
    
    four_sit_result_nn_ITD.append(per_correct_ITD)
    four_sit_result_nn_noITD.append(per_correct_noITD)
        
                 
#print ('SVM:')  
################################################################################################
#print ('C: ' + str(np.mean(four_sit_result[2])) + ' std: ' + str(np.std(four_sit_result[2])))
#print ('A: ' + str(np.mean(four_sit_result[0])) + ' std: ' + str(np.std(four_sit_result[0])))
#print ('D: ' + str(np.mean(four_sit_result[3])) + ' std: ' + str(np.std(four_sit_result[3])))
#print ('B: ' + str(np.mean(four_sit_result[1])) + ' std: ' + str(np.std(four_sit_result[1])))
################################################################################################
#print ('DT:')
################################################################################################
#print ('C: ' + str(np.mean(four_sit_result_DT[2])) + ' std: ' + str(np.std(four_sit_result_DT[2])))
#print ('A: ' + str(np.mean(four_sit_result_DT[0])) + ' std: ' + str(np.std(four_sit_result_DT[0])))
#print ('D: ' + str(np.mean(four_sit_result_DT[3])) + ' std: ' + str(np.std(four_sit_result_DT[3])))
#print ('B: ' + str(np.mean(four_sit_result_DT[1])) + ' std: ' + str(np.std(four_sit_result_DT[1])))
###############################################################################################
print ('nn_ITD:')
###############################################################################################
print ('C: ' + str(np.mean(four_sit_result_nn_ITD[2])))
print ('A: ' + str(np.mean(four_sit_result_nn_ITD[0])))
print ('D: ' + str(np.mean(four_sit_result_nn_ITD[3])))
print ('B: ' + str(np.mean(four_sit_result_nn_ITD[1])))

print ('nn_noITD:')
###############################################################################################
print ('C: ' + str(np.mean(four_sit_result_nn_noITD[2])))
print ('A: ' + str(np.mean(four_sit_result_nn_noITD[0]))) 
print ('D: ' + str(np.mean(four_sit_result_nn_noITD[3])))
print ('B: ' + str(np.mean(four_sit_result_nn_noITD[1])))


## do PCA analysis after wavelet anaylis
# inital 8 PC
# plot top 2 PCA component in 2d graph 

#pca_a = PCA(n_components=8)
#pca_a.fit(wave_coef_matrix_training_feature_selected)
#eigen_a=pca_a.explained_variance_ratio_
#eigen_vector_a=pca_a.components_
#blvg_trace_r_a = pca_a.fit(wave_coef_matrix_training_feature_selected).transform(wave_coef_matrix_training_feature_selected)
#
#x = 0
#y = 5
#plt.plot(blvg_trace_r_a[SSN_index,x], blvg_trace_r_a[SSN_index,y],'go')
#plt.plot(blvg_trace_r_a[speech_index,x], blvg_trace_r_a[speech_index,y],'ro')
