import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# --------------------------------------------------
# 1. Implement AᵀA
# --------------------------------------------------

A = np.array([[1, 2, 3],
              [4, 5, 6]])

AT_A = A.T @ A

print("AᵀA =")
print(AT_A)


# --------------------------------------------------
# 2. y = 2x1 + 3
# --------------------------------------------------

x1 = np.linspace(-100, 100, 100)
y = 2*x1 + 3

plt.figure(figsize=(6,4))
plt.plot(x1, y)
plt.title("y = 2x1 + 3")
plt.xlabel("x1")
plt.ylabel("y")
plt.grid(True)
plt.show()


# --------------------------------------------------
# 3. y = 2x1² + 3x1 + 4
# --------------------------------------------------

x1 = np.linspace(-10, 10, 100)
y = 2*(x1**2) + 3*x1 + 4

plt.figure(figsize=(6,4))
plt.plot(x1, y)
plt.title("y = 2x1² + 3x1 + 4")
plt.xlabel("x1")
plt.ylabel("y")
plt.grid(True)
plt.show()


# --------------------------------------------------
# 4. Gaussian PDF
# mean = 0, sigma = 15
# --------------------------------------------------

mean = 0
sigma = 15

x = np.linspace(-100, 100, 100)
pdf = norm.pdf(x, mean, sigma)

plt.figure(figsize=(6,4))
plt.plot(x, pdf)
plt.title("Gaussian PDF")
plt.xlabel("x")
plt.ylabel("PDF")
plt.grid(True)
plt.show()


# --------------------------------------------------
# 5. y = x1²
# Derivatives at x1 = -5, -3, 0, 3, 5
# --------------------------------------------------

x1 = np.linspace(-10, 10, 100)
y = x1**2

plt.figure(figsize=(6,4))
plt.plot(x1, y)
plt.title("y = x1²")
plt.xlabel("x1")
plt.ylabel("y")
plt.grid(True)
plt.show()

# derivative of x² is 2x
points = np.array([-5, -3, 0, 3, 5])
derivatives = 2 * points

print("\nDerivative values:")
for p, d in zip(points, derivatives):
    print(f"x1 = {p}, dy/dx = {d}")

# value where y = 0
print("\nFunction y = x1² becomes zero at x1 = 0")

print("Inference:")
print("x1 = 0 is the minimum point of the function.")
print("Derivative is negative before 0 and positive after 0.")


# --------------------------------------------------
# 6. y = 2x1 + 3x2 + 3x3 + 4
# Gradient computation
# --------------------------------------------------

def y_function(x1, x2, x3):
    return 2*x1 + 3*x2 + 3*x3 + 4

# Gradient:
# ∂y/∂x1 = 2
# ∂y/∂x2 = 3
# ∂y/∂x3 = 3

print("\nGradient of y:")
print("[2, 3, 3]")

sample_points = [
    (1, 2, 3),
    (0, 1, 1),
    (2, 0, 1)
]

print("\nFunction values at sample points:")
for x1, x2, x3 in sample_points:
    y = y_function(x1, x2, x3)
    print(f"({x1}, {x2}, {x3}) -> y = {y}")


# --------------------------------------------------
# 7. Linear model: Xθ
# --------------------------------------------------

theta = np.array([[2],
                  [3],
                  [3]])

X = np.array([
    [1, 0, 2],
    [0, 1, 1],
    [2, 1, 0],
    [1, 1, 1],
    [0, 2, 1]
])

result = X @ theta

print("\nXθ =")
print(result)