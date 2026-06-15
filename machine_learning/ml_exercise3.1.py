from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import pandas as pd

def load_data():
    # Read ODS file
    data = pd.read_csv("sarrtouris_csv.csv")

    print(data.head())
    print(data.columns)

    # Features (first 3 columns)
    X = data.iloc[:, 2:3]

    # Target (4th column)
    y = data.iloc[:, 3]

    return X, y

def mymain():

    X, y = load_data()

    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=0.20,
        random_state=999
    )

    print("----------TRAINING----------")
    print("N =", len(X))

    model = LinearRegression()
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)

    r2 = r2_score(y_test, y_pred)

    print("R² Score =", round(r2, 4))
    print("Intercept =", model.intercept_)
    print("Coefficients =", model.coef_)

    print("Done")

if __name__ == "__main__":
    mymain()