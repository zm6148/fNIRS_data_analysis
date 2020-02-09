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
# all region answer sheet
filename_path = 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\data\\blvg_trace exp1\\combined_all_answer\\all_answer_sheet.mat'
#filename_path = 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\data\\blvg_trace exp2\\answer\\answer_sheet.mat'
mat_answer = scipy.io.loadmat(filename_path)
answer_sheet_all = mat_answer['all_answer_sheet']

# CA combined answer sheet
filename_path = 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\data\\blvg_trace exp1\\combined_all_answer\\CA_answer_sheet.mat'
#filename_path = 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\data\\blvg_trace exp2\\answer\\answer_sheet.mat'
mat_answer = scipy.io.loadmat(filename_path)
answer_sheet_CA = mat_answer['CA_answer_sheet']

# DB combined answer sheet
filename_path = 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\data\\blvg_trace exp1\\combined_all_answer\\DB_answer_sheet.mat'
#filename_path = 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\data\\blvg_trace exp2\\answer\\answer_sheet.mat'
mat_answer = scipy.io.loadmat(filename_path)
answer_sheet_DB = mat_answer['DB_answer_sheet']

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

## add 2 feaures that corresponds to subject or brain region
# for the subject encoded as 1-of-k
# 14 zeros, first subject change the first 0 to 1 and so on so forwrh

# speech part
#speech_zero_matrix = np.zeros((140, 14))
#starting_index = 0
#for ii in range(0,14): 
#    starting_index = ii * 10
#    speech_zero_matrix[starting_index:starting_index+10, ii] = 1
#   
## SSN part
#SSN_zero_matrix = np.zeros((140, 14))
#starting_index = 0
#for ii in range(0,14):
#    starting_index = ii * 10
#    SSN_zero_matrix[starting_index:starting_index+10, ii] = 1
#    
## combine both part
#subject_matrix = np.concatenate((speech_zero_matrix, SSN_zero_matrix), axis=0)


## if we do the combine of brain regions
# let's set region code in the order as follows
#     C            A          D            B
#     1            0          0            0
# left_tgPCS; right_tgPCS; left_cIFS; right_cIFS 
# and encoded as 1-of-k as well

# combine all regions
all_region_matrix = np.zeros((answer_sheet_all.shape[0], 4))
starting_index = 0
for ii in range(0,4):
    starting_index = ii * 280
    all_region_matrix[starting_index:starting_index+280, ii] = 1
    
# combine CA (tgPCS)
CA_region_matrix = np.zeros((answer_sheet_CA.shape[0], 4))
starting_index = 0
for ii in range(0,2):
    starting_index = ii * 280
    CA_region_matrix[starting_index:starting_index+280, ii] = 1
    
# combine DB (cIFS)
DB_region_matrix = np.zeros((answer_sheet_DB.shape[0], 4))
starting_index = 0
for ii in range(2,4):
    starting_index = (ii-2) * 280
    DB_region_matrix[starting_index:starting_index+280, ii] = 1

    
# path to save all jj of 100 classifilers
jj_path = 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\trained_classifiers\\100_exp1_MNN_with_all_region'

# for jj times repeat train 100 classifier on each recording site.
for jj in range(0,1):
    
    # for this jj save trained 100 clasifier at
    jj_path_each = jj_path + '\\' + str(jj)
    
    # repeat this whole process k times
    # to aclculate the averrage accuracy

    four_sit_result = []
    four_sit_result_DT = []
    four_sit_result_nn = []
    four_sit_loss = []
    for filename in os.listdir(dirname):
        
        
        # load matlab file and save filename
        print(filename)
        filename_list.append(filename)
        site_name = filename.split('.')[0]
        
        filename_path=dirname+'/'+filename
        mat = scipy.io.loadmat(filename_path)
        data_matrix = mat[site_name]
        
        # based on file name decide which answer sheet and region feature to use
        if filename == 'all_matrix.mat':
            answer_sheet = answer_sheet_all;
            region_matrix = all_region_matrix      
        elif filename == 'CA_matrix.mat':
            answer_sheet = answer_sheet_CA;
            region_matrix = CA_region_matrix  
            # Do the other thing
        elif filename == 'DB_matrix.mat':
            # Do yet another thing
            answer_sheet = answer_sheet_DB;
            region_matrix = DB_region_matrix  

        
        # divide to 2 parts: training and testing
        # random divide for n times
        # each time to training and classification
        all_per_correct = []
        all_per_correct_DT = []
        all_per_correct_nn = []
        all_loss = []
        repeat_num = 100
        for n in range(0,repeat_num):
            
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
            
            ## build speech_zero_matrix that is training selected
            #subject_matrix_training = subject_matrix[training_index,:]
            #subject_matrix_testig = subject_matrix[testing_index,:]
            
            ## build region matreix that is training selected
            region_matrix_training = region_matrix[training_index,:]
            region_matrix_testig = region_matrix[testing_index,:]
            
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
        
            ## using 25 features extracted
            wave_coef_matrix_training_feature_selected = wave_coef_matrix_training[:,selected_feature_index]
            
            ## combibe wave_coef_matrix_training_feature_selected with subject columns
            wave_coef_matrix_training_feature_selected_regions = np.concatenate((wave_coef_matrix_training_feature_selected, region_matrix_training), axis=1)
            
            X = wave_coef_matrix_training_feature_selected_regions
            y = training_answer
            
            ## train svm
            #clf = svm.SVC()
            ## X is the training data set
            ## y is the class label
            #X = wave_coef_matrix_training_feature_selected_subject_columns
            #y = training_answer
            #clf.fit(X, y)  
            #
            ## train DT
            #clf_DT = tree.DecisionTreeClassifier()
            #clf_DT = clf_DT.fit(X, y)
            #
            ## use clf to perdict 
            ## use the same feature as training data
            #out_put = clf.predict(wave_coef_matrix_testig[:,selected_feature_index])
            #out_put_DT = clf_DT.predict(wave_coef_matrix_testig[:,selected_feature_index])
            #
            ## see how good it is
            #per_correct = fun.predict_accuracy(testing_answer, out_put)
            #per_correct_DT = fun.predict_accuracy(testing_answer, out_put_DT)
            ## add the results to all_per_correct
            #all_per_correct.append(per_correct)
            #all_per_correct_DT.append(per_correct_DT)
            
            ## NN classifier
            bestseed = 0
            bestloss = 10000000
            for i in range(0, 10, 1):
                r = random.randint(0,10000000)    
                clf_nn = MLPClassifier(solver='lbfgs', activation='logistic', hidden_layer_sizes=(20,10,5),random_state=r,max_iter=1000)
                clf_nn.fit(X,y)
                
                if(clf_nn.loss_ < bestloss):
                    bestloss = clf_nn.loss_
                    bestseed = r
                
            print(clf_nn.loss_, bestloss)
            
            # save loss
            all_loss.append(bestloss)
            
            clf_nn = MLPClassifier(solver='lbfgs', activation='logistic', hidden_layer_sizes=(20,10,5),random_state=bestseed,max_iter=1000)
            clf_nn.fit(X,y)
            
            # before making perdictions
            # add the subject column to testing matrix as well
            wave_coef_matrix_testig_subject = np.concatenate((wave_coef_matrix_testig[:,selected_feature_index], region_matrix_testig), axis=1)
            prediction = clf_nn.predict(wave_coef_matrix_testig_subject)
            per_correct_nn = fun.predict_accuracy(testing_answer, prediction)
            all_per_correct_nn.append(per_correct_nn)
            
            ## each time create a folder and save trained classifer in it.
            path  = jj_path_each + '\\' + site_name + '\\'+ str(n);
            try: 
                os.makedirs(path)
            except OSError:
                if not os.path.isdir(path):
                    raise
                    
            #joblib.dump(clf, path+'\\clf_svm.pkl') 
            #joblib.dump(clf_DT, path+'\\clf_DT.pkl') 
            joblib.dump(clf_nn, path+'\\clf_nn.pkl') 
            
            # also slave selected feature index
            np.save(path+'\\feature_index.npy', selected_feature_index)
            
        
        all_per_correct = np.asarray(all_per_correct)
        all_per_correct_DT = np.asarray(all_per_correct_DT)
        all_per_correct_nn = np.asarray(all_per_correct_nn)
        all_loss = np.asarray(all_loss)
        
        four_sit_result.append(all_per_correct)
        four_sit_result_DT.append(all_per_correct_DT)
        four_sit_result_nn.append(all_per_correct_nn)
        four_sit_loss.append(all_loss)
        
    # save the accuracy for later
    np.save(jj_path_each +'\\accuracy.npy', four_sit_result_nn)
    # save the best lose for later
    np.save(jj_path_each +'\\loss.npy', four_sit_loss)

    
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
print ('nn:')
###############################################################################################
print ('All combined: ' + str(np.mean(four_sit_result_nn[0])) + ' std: ' + str(np.std(four_sit_result_nn[0])))
print ('CA: ' + str(np.mean(four_sit_result_nn[1])) + ' std: ' + str(np.std(four_sit_result_nn[1])))
print ('DB: ' + str(np.mean(four_sit_result_nn[2])) + ' std: ' + str(np.std(four_sit_result_nn[2])))
#print ('B: ' + str(np.mean(four_sit_result_nn[1])) + ' std: ' + str(np.std(four_sit_result_nn[1])))

### save trained classifier for later use on exp2
## save trained classifier to a file
#joblib.dump(clf, 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\clf_svm.pkl') 
#joblib.dump(clf_DT, 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\clf_DT.pkl') 
#joblib.dump(clf_nn, 'C:\\Users\\mz86\\Desktop\\fNIRS ML\\clf_nn.pkl') 
#
## also slave selected feature index
#np.save('C:\\Users\\mz86\\Desktop\\fNIRS ML\\feature_index.npy', selected_feature_index)

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
