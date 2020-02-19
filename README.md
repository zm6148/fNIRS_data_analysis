# fNIRS_data_analysis


# fNIRS_ML.ipynb is the colab file with ML on fnirs data procedures
https://colab.research.google.com/drive/1GAEBau--6o1TsjpCUo6DFkH9camo3syn#scrollTo=BukfmkkTuWbU

# data preprocessing
 matlab code for fNIRS analysis
 Each experiment is located in separate folder.
 
 data used can be found on
 https://drive.google.com/drive/folders/1v6MDBZKNJeczF59sMhIQATFdXPD2QVkx?usp=sharing
 
 
Those scrips were written in MatLab R2016a. And “Homer2” is the folder for Homer2, a set of matlab scripts used for analyzing fNIRS data (http://homer-fnirs.org/). Add this folder to matlab path before running analysis scripts. Also matlab wavelet toolbox is required for wavelet analysis used in the scripts.

“Subject keypress data” folder contains the subject key presses during the auditory attention task and passive listening task for both experiments 1 and subject key presses during the auditory attention task in experiment 2 stored in folders named: “experiment 1” and “experiment 2” respectively.

1)	“calculate_ percentage_correct.m” is the matlab analysis code to calculate the percentage correct of subject response. For example:
Load the subject response file “TIC_2_conditions_all_key_press_block.mat” and run “calculate_ percentage_correct.m” will output the accuracy for subject TIC.

2)	“calculate_percentage_correct_simulation.m” is the script to simulate the percentage correct if subject responded randomly. The average percentage correct after 1000000 trails is 16.87%.

“Experiment 1” folder contains all the data and analysis for experiment 1.

1)	“Raw Recordings” is the folder for original finrs recordings for experiment 1. It has two sub-folders one for each of the two task conditions: Speech or SSN. Those were in the .nirs format.

2)	Those original recordings were change to .mat format and placed in “…\Fnirs data and analysis\Experiment 1\Data\Speech\Unprocessed”  and “…\Fnirs data and analysis\Experiment 1\Data\SSN\Unprocessed”

3)	Run the .m file named “data_preprocess_exp1.m” to begin the preprocessing steps for all recordings. Recordings from Speech condition and SSN condition can be loaded separately by choosing different directory at the beginning of this .m file. 

4)	Once the preprocessing is done, the processed data will be stored in “…\Fnirs data and analysis\Experiment 1\Data\Speech\Processed” and “…\Fnirs data and analysis\Experiment 1\Data\SSN\Processed”

5)	To being GLM analysis run the .m file named: “data_analysis_GLM_exp1.m”. 2 figures will appear, figure 1 is the GLM beta value for the breath holding task at each recording site:
C (left tgPCS);  A (right tgPCS); D (left cIFS);   (right cIFS). Figure 2 is the GLM beta value for the perceptual under 2 conditons: Speech and SSN. The GLM results will be stored in “GLM_data_exp1.mat”  file. 

6)	Too see the block average trace of HbO, run the .m file named “data_analysis_HbO_trace_exp1.m” and choose at the beginning of code to set block average from which condition to see.


“Experiment 2” folder contains all the data and analysis for experiment 2. It follows a similar  directory structure as “Experiment 1”

1)	“Raw Recordings” is the folder for original finrs recordings for experiment 2. It has two sub-folders one for each of the two task conditions: Speech or SSN. Those were in the .nirs format.

2)	Those original recordings were change to .mat format and placed in “…\Fnirs data and analysis\Experiment 1\ Data\Unprocessed”.

3)	Run the .m file named “data_preprocess_exp2.m” to begin the preprocessing steps for all recordings. Recordings from Speech condition and SSN condition can be loaded separately by choosing different directory at the beginning of this .m file. 

5)	Once the preprocessing is done, the processed data will be stored in “…\Fnirs data and analysis\Experiment 1\Data\ Processed” 

6)	To being GLM analysis run the .m file named: “data_analysis_GLM_exp2.m”. 2 figures will appear, figure 1 is the GLM beta value for the breath holding task at each recording site:
C (left tgPCS);  A (right tgPCS); D (left cIFS);   (right cIFS). Figure 2 is the GLM beta value for the perceptual under 2 conditons: pitch and spatial cue and pitch cue only. The GLM results will be stored in “GLM_data_exp2.mat”  file. 

7)	Too see the block average trace of HbO, run the .m file named “data_analysis_HbO_trace_exp2.m”.  

Example results from Exp 1
1. block average
![GitHub Logo](https://lh3.googleusercontent.com/mC4m-J3bE8xoLZLf-VGacvzjBuoSnf_Vu8RiF9YoPV8d7QlChoeFcZdwok9mczSTwuUL3hfIFdehjUswl1_33gSs8PjjTz6-IJ8APLUpLxLeDMCRnbwfVSs5A4_oKQ4a-g1ezDAG5g=w2400)
2. GLM
![GitHub Logo](https://lh3.googleusercontent.com/2AtLXjpsQKFcV6JoeIuszot7RBCqbT9cXd2osBZX41h8efZqHpRUVKywvXDzujNYADulz5XD-VL2n8tCx1N-t1u9IWoj74CZQOAlCXlZLFFOP3LWzamTOi3s7Yid-FT3t4t6A2oy1g=w2400)


ML on fnirs data is the folder for ML using fNIRS block average data
1. how to generate data set from limited subjects
![GitHub Logo](https://lh3.googleusercontent.com/dcT8XDDT8OFxplPU3BR8FNh0ZFOp7TayGYAV7bVV8vQ_DSgR9Ppq3lZPnPUXfH1FTP4M2NT_VjPlGw6bjjgfLHa9f3q2-_nbaF2Ori83QFBUaOjXkX1yXZ0_2ZOffyKDwd7PtAbxzQ=w2400)
2. SVM results
![GitHub Logo](https://lh3.googleusercontent.com/IptTuhyG3MwcT0YgIbMZ5z3shnEQe18jPADROBCxOxWWvNLdvcgngVVTPN3EwtUOeDuUiVQw1VaSpwhvifiQ7AeIB7RECrPv9R-28KAFPdq3r5hQQswRAieNW35R_l1csPqlotbOXQ=w2400)

