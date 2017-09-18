import os
from zipfile import ZipFile


files = ['h1_classifiers.py',
         'logistic_regression.py',
         'model_stats.py',
         'softmax.py',
         'h1_util.py']

results = ['logistic_regression_best_result.csv',
           'logistic_regression_two_class_parameter_plot_1_128.png',
           'logistic_regression_validation_accuracy.png',
           'logistic_regression_validation_scores.csv',
           'results_2_vs_7.csv',
           'softmax_best_result.csv',
           'softmax_confusion_matrix.csv',
           'softmax_parameter_plot_1_128.png',
           'softmax_stats_accuracy.csv',
           'softmax_validation_accuracy.png',
           'softmax_validation_scores.csv']

with ZipFile('handin1_upload_files.zip', 'w') as myzip:
    for filename in files:
        myzip.write(filename)
    for filename in results:
        rf = os.path.join('results', filename)
        myzip.write(rf)

