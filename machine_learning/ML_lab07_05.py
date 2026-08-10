import numpy as np
import pandas as pd
from sklearn.datasets import load_diabetes, load_iris
from sklearn.model_selection import train_test_split
from xgboost import XGBRegressor, XGBClassifier
from sklearn.metrics import mean_squared_error, accuracy_score, classification_report, confusion_matrix
# Load the Diabetes dataset
diabetes = load_diabetes()
X = diabetes.data
y = diabetes.target

# Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2,
random_state=42)

# Initialize the XGBoost Regressor
xgb_regressor = XGBRegressor(objective='reg:squarederror',
n_estimators=100, learning_rate=0.1, max_depth=3, random_state=42)

# Train the model on the training data
xgb_regressor.fit(X_train, y_train)

# Predict using the trained model on the testing set
y_pred = xgb_regressor.predict(X_test)

# Evaluate the model performance
mse = mean_squared_error(y_test, y_pred)
print(f"Mean Squared Error: {mse:.4f}")

# Load the Iris dataset
iris = load_iris()
X = iris.data
y = iris.target

# Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2,
random_state=42)

# Initialize the XGBoost Classifier
xgb_classifier = XGBClassifier(objective='multi:softmax', num_class=3,
n_estimators=100, learning_rate=0.1, max_depth=3, random_state=42)

# Train the model on the training data
xgb_classifier.fit(X_train, y_train)

# Predict using the trained model on the testing set
y_pred = xgb_classifier.predict(X_test)

# Evaluate the model performance
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy: {accuracy:.4f}")

report = classification_report(y_test, y_pred)
conf_matrix = confusion_matrix(y_test, y_pred)

print("Classification Report:")
print(report)

print("Confusion Matrix:")
print(conf_matrix)

