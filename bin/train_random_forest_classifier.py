#!/usr/bin/python

# BUILD AND TRAIN A RANDOM FOREST CLASSIFIER

# IMPORT NECESSARY LIBRARIES
import sys
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report

# LOAD THE DATASET
data_file = 'path/to/your/data.csv'
data = pd.read_csv(data_file)

# SEPARATE FEATURES AND TARGET
x = data.iloc[:, :-1] # Select all columns except the last one
y = data.iloc[:, -1] # Select the last column

# SPLIT THE DATASET INTO TRAINING AND TEST SETS
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)

# INITIALIZE THE RANDOM FOREST CLASSIFIER
rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)

# TRAIN THE CLASSIFIER
rf_classifier.fit(x_train, y_train)

# MAKE PREDICTIONS ON TEST SET
y_pred = rf_classifier.predict(x_test)

# EVALUATE THE CLASSIFIER
accuracy = accuracy_score(y_test, y_pred)
print(f'Accuracy: {accuracy:.2f}')

# PRINT A DETAILED CLASSIFICATION REPORT
report = classification_report(y_test, y_pred)
print('Classification Report:')
print(report)

# OPTIONAL: SAVE THE TRAINED MODEL FOR LATER USE
import joblib
model_file = 'cnv_rf_classifier.pkl'
joblib.dump(rf_classifier, model_file)
print(f'Model saved to {model_file}')