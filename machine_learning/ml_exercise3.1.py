from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import pandas as pd
def load_data():
    data=pd.read_csv('simulated_data_multiple_linear_regression_for_ML.csv')
    X=data.iloc[:,:5]
    y=data.iloc[:,6]
    return X,y
def mymain():

    # Load dataset
    X, y = load_data()

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X,y,test_size=0.20,random_state=999)

    print('----------TRAINING----------')
    print("N=%d" % len(X))

    # Create model
    model = LinearRegression()

    # Train model
    model.fit(X_train, y_train)

    # Prediction
    y_pred = model.predict(X_test)

    # Accuracy
    r2 = r2_score(y_test, y_pred)

    print("r2 score is %0.2f (closer to 1 is good)" % r2)

    print("done")


if __name__ == "__main__":
    mymain()