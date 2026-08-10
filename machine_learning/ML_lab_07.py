from sklearn.tree import DecisionTreeRegressor
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import pandas as pd
import matplotlib.pyplot as plt

# Generate synthetic dataset
X, y = make_regression(
    n_samples=1000,
    n_features=20,
    noise=0.1,
    random_state=42
)

# Split dataset
X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size=0.2,
    random_state=42
)

# Create model
regressor = DecisionTreeRegressor(max_depth=5, random_state=42)

# Train
regressor.fit(X_train, y_train)

# Predict
y_pred = regressor.predict(X_test)

# Evaluate
mse = mean_squared_error(y_test, y_pred)
print(f"Mean Squared Error: {mse:.4f}")

# Feature Importance
plt.figure(figsize=(10,6))

feature_importances = pd.Series(
    regressor.feature_importances_,
    index=[f"Feature {i}" for i in range(X.shape[1])]
)

feature_importances.nlargest(10).sort_values().plot(kind='barh')

plt.title("Top 10 Feature Importances")
plt.xlabel("Importance")
plt.tight_layout()
plt.show()

