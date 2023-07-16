import numpy as np
from field import Field17


class Polynomial:
    """A naive implementation of polynomials in coefficient form.

    The i-th term is the coefficient of the i-th power, eg.
    f(x) = 2x^2 + 4x + 1 will be represented as [1, 4, 2].
    """
    coefficients = []

    def __init__(self, coefficients):
        self.coefficients = [Field17(coefficient) for coefficient in coefficients]

    def __len__(self):
        return len(self.coefficients)

    def degree(self):
        if len(self.coefficients) == 0:
            return -1

        degree = 0
        for i in range(len(self.coefficients)):
            if self.coefficients[i] != 0:
                degree = i

        return degree

    def leading_coefficient(self):
        return self.coefficients[self.degree()]

    def __repr__(self) -> str:
        return f"Polynomial([{', '.join(str(c) + 'x^' + str(i) for i, c in enumerate(self.coefficients))}])"

    def __radd__(self, other):
        return self.__add__(other)

    def __add__(self, other):
        if isinstance(other, int) or isinstance(other, Field17):
            coefficients = [coeff for coeff in self.coefficients]
            coefficients[0] = coefficients[0] + other
            return Polynomial(coefficients)
        if self.degree() == -1:
            return other
        elif other.degree() == -1:
            return self
        coeffs = [0] * max(len(self.coefficients), len(other.coefficients))
        for i in range(len(self.coefficients)):
            coeffs[i] = coeffs[i] + self.coefficients[i]
        for i in range(len(other.coefficients)):
            coeffs[i] = coeffs[i] + other.coefficients[i]
        return Polynomial(coeffs)

    def __neg__(self):
        return Polynomial([-c for c in self.coefficients])

    def __sub__(self, other):
        return self.__add__(-other)

    def __getitem__(self, i):
        return self.coefficients[i]

    def __iter__(self):
        return iter(self.coefficients)

    def __eq__(self, other):
        if self.degree() != other.degree():
            return False
        if self.degree() == -1:
            return True
        return self.coefficients == other.coefficients

    def __rmul__(self, other):
        return self.__mul__(other)

    def __mul__(self, other):
        if isinstance(other, Field17):
            return Polynomial(
                [coefficient * other for coefficient in self.coefficients]
            )

        if isinstance(other, int):
            return Polynomial(
                [coefficient * other for coefficient in self.coefficients]
            )
        degree = len(self.coefficients) + len(other.coefficients) - 1

        new_coeffs = [0 for _ in range(degree)]
        for i, coeff in enumerate(self.coefficients):
            for j, coeff_other in enumerate(other.coefficients):
                new_coeffs[i + j] += coeff * coeff_other

        while new_coeffs[-1] == 0:
            new_coeffs.pop()
        return Polynomial([coeff for coeff in new_coeffs])

    def __truediv__(self, other):
        quo, rem = Polynomial.divide(self, other)

        assert (
            rem.degree() == -1
        ), "cannot perform polynomial division because remainder is not zero"
        return quo

    def divide(numerator, denominator):
        if len(denominator) == 1:
            return None

        if len(numerator) < len(denominator):
            return (Polynomial([]), numerator)

        remainder = Polynomial([n for n in numerator.coefficients])

        quotient_coefficients = [
            0 for _ in range(numerator.degree() - denominator.degree() + 1)
        ]

        for i in range(numerator.degree() - denominator.degree() + 1):
            if remainder.degree() < denominator.degree():
                break
            coefficient = (
                remainder.leading_coefficient() / denominator.leading_coefficient()
            )
            shift = remainder.degree() - denominator.degree()
            subtractee = Polynomial([0] * shift + [coefficient]) * denominator
            quotient_coefficients[shift] = coefficient
            remainder = remainder - subtractee
        quotient = Polynomial(quotient_coefficients)
        return quotient, remainder

    def synthetic_divide(self, divisor):
        """Simpler way than long division to divide polynomials."""
        assert (
            divisor.degree() == 1
        ), f"Divisor polynomial must be of degree 1. Received: {divisor}"

        divisor = divisor.coefficients[0]
        for i in range(1, len(self.coefficients) - 1):
            i = self.degree() - i
            result = self.coefficients[i] + self.coefficients[i + 1] * divisor
            self.coefficients[i] = result

        return Polynomial(self.coefficients[1:])

    def synthetic_divide_similar(self, other):
        """Divides 2 polynomials using the method similar to synthetic division in the blog.

        Work from right to left (highest power to lowest power), take each coefficient
        and add the coefficient directly below. If no coefficient, take it as 0.
        """
        new_coefficients = []

        for i, coefficient in enumerate(self[::-1]):
            new_coefficient = coefficient
            if i > other.degree() - 1:
                new_coefficient = coefficient + new_coefficients[i - other.degree()]
            new_coefficients.append(new_coefficient)

        return Polynomial(new_coefficients[::-1][other.degree() : :])

    def __mod__(self, other):
        quo, rem = Polynomial.divide(self, other)
        return rem

    def eval(self, x):
        summed = 0
        for i, coefficient in enumerate(self.coefficients):
            summed += coefficient * (pow(x, i))

        return summed


def interpolate_poly_degree_3(xs, B):
    """Interpolates a degree-3 polynomial using a matrix equation.
    Given f(x) = d + cx + bx^2 + ax^3, express the coefficients and x values
    as matrices:

    Ax = B => x = A^-1 * B

    where A is our roots of unity and B is our desired outputs.

    Taking the inverse of A allows us to solve the system of equations.
    """

    A = np.array([[x**d for d in range(4)] for x in xs])

    A_inverse = inverse_matrix(A)

    return Polynomial([Field17(i) for i in A_inverse @ B])


def inverse_matrix(matrix):
    n = matrix.shape[0]
    augmented_matrix = np.hstack((matrix, np.identity(n, dtype=int)), dtype=int)

    for i in range(n):
        # Find the modular multiplicative inverse of the diagonal element
        diagonal_element = augmented_matrix[i, i]
        inverse_diagonal = pow(
            int(diagonal_element), -1, 17
        )  # Modular multiplicative inverse

        # Multiply the entire row by the modular multiplicative inverse
        augmented_matrix[i] = augmented_matrix[i] * inverse_diagonal

        for j in range(n):
            if j != i:
                # Subtract a multiple of the current row to make elements below the diagonal zero
                factor = augmented_matrix[j, i]
                augmented_matrix[j] = augmented_matrix[j] - factor * augmented_matrix[i]

    # Extract the inverse matrix from the augmented matrix
    inverse = augmented_matrix[:, n:]

    return inverse.astype(int)
