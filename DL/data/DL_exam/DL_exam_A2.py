import numpy as np
import matplotlib.pyplot as plt


z = np.linspace(-10, 10, 100)

# Sigmoid function and its derivative
def sigmoid(z):
    return 1 / (1 + np.exp(-z))

def sigmoid_derivative(z):
    s = sigmoid(z)
    ds = s * (1 - s)
    return ds

# Tanh function and its derivative
def tanh(z):
    return np.tanh(z)

def tanh_derivative(z):
    return 1 - np.tanh(z) ** 2

# ReLU function and its derivative
def relu(z):
    return np.maximum(0, z)

def relu_derivative(z):
    return np.where(z > 0, 1, 0)

def main():

    plt.plot(z, sigmoid(z), label='Sigmoid')
    plt.plot(z, sigmoid_derivative(z), label='Sigmoid Derivative')
    plt.legend()
    plt.plot(z, tanh(z), label='Tanh')
    plt.plot(z, tanh_derivative(z),label='Tanh Derivative')
    plt.legend()
    plt.plot(z, relu(z), label='ReLU')
    plt.plot(z, relu_derivative(z), label='ReLU Derivative')
    plt.title('All activation function  and its Derivative')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    # print("this is the beginning of my program")
    main()

