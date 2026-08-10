from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt
import pandas as pd

# Load the Iris dataset
iris = load_iris()
X = iris.data
y = iris.target

# Split the dataset
X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size=0.2,
    random_state=42
)

# Create Decision Tree Classifier
classifier = DecisionTreeClassifier(max_depth=5, random_state=42)

# Train the model
classifier.fit(X_train, y_train)

# Predict
y_pred = classifier.predict(X_test)

# Evaluate
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy: {accuracy:.4f}")

print("\nClassification Report:")
print(classification_report(y_test, y_pred))

print("\nConfusion Matrix:")
print(confusion_matrix(y_test, y_pred))

# Plot Feature Importance
feature_importances = pd.Series(
    classifier.feature_importances_,
    index=iris.feature_names
)

plt.figure(figsize=(8,5))

feature_importances.sort_values().plot(kind='barh')

plt.title("Feature Importances")
plt.xlabel("Importance")
plt.tight_layout()
plt.show()