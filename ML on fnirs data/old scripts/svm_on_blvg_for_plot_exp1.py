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
filename_path = 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\data\\blvg_trace exp1\\answer\\answer_sheet.mat'
#filename_path = 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\data\\blvg_trace exp2\\answer\\answer_sheet.mat'
mat_answer = scipy.io.loadmat(filename_path)
answer_sheet = mat_answer['answer_sheet']

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
four_sit_result_nn = []
for filename in os.listdir(dirname):
    
    # load matlab file and save filename
    print(filename)
    filename_list.append(filename)
    site_name = filename.split('.')[0]
    
    filename_path=dirname+'/'+filename
    mat = scipy.io.loadmat(filename_path)
    data_matrix = mat[site_name]
    
    # divide to 2 parts: training and testing
    # random divide for n times
    # each time to training and classification
    all_per_correct = []
    all_per_correct_DT = []
    all_per_correct_nn = []
    repeat_num = 100
    for ii in range(0,repeat_num):
        
        row_num = np.shape(data_matrix)[0]
        all_index = np.array(range(0, row_num))
        
        # generate training index and testing index
        training_index = random.sample(all_index, int(row_num * 0.9))
        testing_index = np.delete(all_index, training_index)
        
        ## select training set and testing set
        #training_set = data_matrix[training_index]
        #testing_set = data_matrix[testing_index]
        
        # create training and testing answer
        testing_answer = answer_sheet[testing_index]
        training_answer = answer_sheet[training_index]
        
        speech_index = np.where(training_answer == 1)[0]
        SSN_index = np.where(training_answer == 0)[0]
        
        
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
        
        ## to train on training set only
        wave_coef_matrix_training = wave_coef_matrix[training_index,:]
        wave_coef_matrix_testig = wave_coef_matrix[testing_index,:]
        
        ## from wave_coef_matrix choose 10 columns least likely to be random
        # using Kolmogorov-Smirnov (KS) test
        # http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.653.6292&rep=rep1&type=pdf
        all_feature_test_result = []
        for ii in range(0, wave_coef_matrix_training.shape[1]):
            feature = wave_coef_matrix_training[:,ii].astype(float)
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
    
        ## do PCA analysis after wavelet anaylis
        ## inital 8 PC
        ## plot top 2 PCA component in 2d graph 
        #pca_a = PCA(n_components=8)
        #pca_a.fit(wave_coef_matrix)
        #eigen_a=pca_a.explained_variance_ratio_
        #eigen_vector_a=pca_a.components_
        #blvg_trace_r_a = pca_a.fit(wave_coef_matrix).transform(wave_coef_matrix)
    
        ## start training SVM using 25 features extracted
        wave_coef_matrix_training_feature_selected = wave_coef_matrix_training[:,selected_feature_index]
        
        # train svm
        clf = svm.SVC()
        # X is the training data set
        # y is the class label
        X = wave_coef_matrix_training_feature_selected
        y = training_answer
        clf.fit(X, y)  
        
        # train DT
        clf_DT = tree.DecisionTreeClassifier()
        clf_DT = clf_DT.fit(X, y)
        
        # use clf to perdict 
        # use the same feature as training data
        out_put = clf.predict(wave_coef_matrix_testig[:,selected_feature_index])
        out_put_DT = clf_DT.predict(wave_coef_matrix_testig[:,selected_feature_index])
        
        # see how good it is
        per_correct = fun.predict_accuracy(testing_answer, out_put)
        per_correct_DT = fun.predict_accuracy(testing_answer, out_put_DT)
        # add the results to all_per_correct
        all_per_correct.append(per_correct)
        all_per_correct_DT.append(per_correct_DT)
        
        ## NN classifier
        bestseed = 0
        bestloss = 10000000
        for i in range(0, 10, 1):
            r = random.randint(0,10000000)    
            clf_nn = MLPClassifier(solver='lbfgs', activation='logistic', hidden_layer_sizes=(15),random_state=r,max_iter=1000)
            clf_nn.fit(X,y)
            
        if(clf_nn.loss_ < bestloss):
            bestloss = clf_nn.loss_
            bestseed = r
        print(clf_nn.loss_, bestloss)
        
        clf_nn = MLPClassifier(solver='lbfgs', activation='logistic', hidden_layer_sizes=(15),random_state=bestseed,max_iter=1000)
        clf_nn.fit(X,y)
        
        prediction = clf_nn.predict(wave_coef_matrix_testig[:,selected_feature_index])
        per_correct_nn = fun.predict_accuracy(testing_answer, prediction)
        all_per_correct_nn.append(per_correct_nn)
        
        
    
    all_per_correct = np.asarray(all_per_correct)
    all_per_correct_DT = np.asarray(all_per_correct_DT)
    all_per_correct_nn = np.asarray(all_per_correct_nn)
    
    four_sit_result.append(all_per_correct)
    four_sit_result_DT.append(all_per_correct_DT)
    four_sit_result_nn.append(all_per_correct_nn)
    
print ('SVM:')  
###############################################################################################
print ('C: ' + str(np.mean(four_sit_result[2])) + ' std: ' + str(np.std(four_sit_result[2])))
print ('A: ' + str(np.mean(four_sit_result[0])) + ' std: ' + str(np.std(four_sit_result[0])))
print ('D: ' + str(np.mean(four_sit_result[3])) + ' std: ' + str(np.std(four_sit_result[3])))
print ('B: ' + str(np.mean(four_sit_result[1])) + ' std: ' + str(np.std(four_sit_result[1])))
###############################################################################################
print ('DT:')
###############################################################################################
print ('C: ' + str(np.mean(four_sit_result_DT[2])) + ' std: ' + str(np.std(four_sit_result_DT[2])))
print ('A: ' + str(np.mean(four_sit_result_DT[0])) + ' std: ' + str(np.std(four_sit_result_DT[0])))
print ('D: ' + str(np.mean(four_sit_result_DT[3])) + ' std: ' + str(np.std(four_sit_result_DT[3])))
print ('B: ' + str(np.mean(four_sit_result_DT[1])) + ' std: ' + str(np.std(four_sit_result_DT[1])))
###############################################################################################
print ('nn:')
###############################################################################################
print ('C: ' + str(np.mean(four_sit_result_nn[2])) + ' std: ' + str(np.std(four_sit_result_nn[2])))
print ('A: ' + str(np.mean(four_sit_result_nn[0])) + ' std: ' + str(np.std(four_sit_result_nn[0])))
print ('D: ' + str(np.mean(four_sit_result_nn[3])) + ' std: ' + str(np.std(four_sit_result_nn[3])))
print ('B: ' + str(np.mean(four_sit_result_nn[1])) + ' std: ' + str(np.std(four_sit_result_nn[1])))

## save trained classifier for later use on exp2
# save trained classifier to a file
joblib.dump(clf, 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\clf_svm.pkl') 
joblib.dump(clf_DT, 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\clf_DT.pkl') 
joblib.dump(clf_nn, 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\clf_nn.pkl') 

# also slave selected feature index
np.save('C:\\Users\\mz86\\Desktop\\fNIRS ML\\feature_index.npy', selected_feature_index)

# and later you can load it
#clf = joblib.load('filename.pkl')




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
