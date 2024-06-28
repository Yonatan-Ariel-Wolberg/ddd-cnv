#!/usr/bin/python

# LOAD THE NECESSARY LIBRARIES
import pandas as pd
from sklearn.metrics import accuracy_score, classification_report
import joblib

# LOAD THE NEW CNV DATA INTO PANDAS DATAFRAME
new_data = pd.read_csv('new_cnv_data.csv')

# ENSURE THE NEW DATA HAS THE SAME FEATURES AS THE TRAINING DATA
# Assuming your training data had a target column 'label' which we drop for predictions
x_new = new_data.drop('label', axis=1, errors='ignore') # Drop the label column if it exists

# LOAD THE TRAINED MODEL
loaded_model = joblib.load('cnv_rf_classifier.pkl')

# MAKE PREDICTIONS ON THE NEW DATA
new_predictions = loaded_model.predict(x_new)

# IF THE NEW CNV DATA HAS GROUND TRUTH LABELS, EVALUATE THE PREDICTIONS
if 'label' in new_data.columns:
    y_new = new_data['label']
    print("Classification Report:")
    print(classification_report(y_new, new_predictions))

print("Accuracy Score:")
print(accuracy_score(y_new, new_predictions))

# SAVE PREDICTIONS TO NEW CSV FILE
output = new_data.copy()
output['predictions'] = new_predictions
output.to_csv('validated_cnv_data.csv', index=False)

print("Predictions saved to 'validated_cnv_data.csv'")
