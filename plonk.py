import unittest
from enum import Enum


# Elliptic curve: y^2 = x^3 + 3


class Field17:
    modulus = 17

    def __init__(self, value):
        self.value = value % self.modulus


class PrimeField:
    modulus = 101

    def __init__(self, value, modulus=None):
        if modulus:
            self.modulus = modulus
        self.value = value % self.modulus

    def __add__(self, other):
        return PrimeField(self.value + other.value)

    def __sub__(self, other):
        return PrimeField(self.value - other.value)

    def __rmul__(self, other):
        return self * other

    def __mul__(self, other):
        if isinstance(other, PrimeField):
            value = other.value
        elif isinstance(other, int):
            value = other
        return PrimeField(self.value * value)

    def __pow__(self, other):
        if isinstance(other, PrimeField):
            value = other.value
        elif isinstance(other, int):
            value = other
        return PrimeField(self.value**value)

    def __truediv__(self, other):
        return self.__div__(other)

    def __div__(self, other):
        if isinstance(other, PrimeField):
            value = other
        elif isinstance(other, int):
            value = PrimeField(other)

        return PrimeField(self.value * value.inv())

    def __eq__(self, other):
        if isinstance(other, PrimeField):
            value = other.value
        elif isinstance(other, int):
            value = other

        return self.value == value

    def inv(self):
        """
        Extended euclidean algorithm to find modular inverses for integers
        """
        # To address a == n edge case.
        # https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-09#section-4
        # inv0(x): This function returns the multiplicative inverse of x in
        # F, extended to all of F by fixing inv0(0) == 0.
        self.value %= self.modulus

        if self.value == 0:
            return 0
        lm, hm = 1, 0
        low, high = self.value % self.modulus, self.modulus
        while low > 1:
            r = high // low
            nm, new = hm - lm * r, high - low * r
            lm, low, hm, high = nm, new, lm, low

        return lm % self.modulus


class Generator:
    def __init__(self, x, y, modulus):
        self.x = x % modulus
        self.y = y % modulus


class StructuredReferenceString:
    values = []

    def __init__(self):
        powers_of_s = [2, 4, 8, 16, 15, 13]

        g = G1(1, 2)
        first = [G1(1, 2)]
        for _ in range(len(powers_of_s)):
            first.append(g.double())

        g2 = G2(36, 31)
        second = [G2(36, 31)]
        doubled = g2.double()
        second.append(doubled)

        self.values = first + second

    def __iter__(self):
        return iter(self.values)

    def __eq__(self, other):
        return all(
            a.x.value == b.x.value and a.y.value == b.y.value
            for (a, b) in zip(self, other)
        )


class G2:
    prime = 101

    def __init__(self, x, y):
        self.x = PrimeField(x)
        self.y = PrimeField(y)

    def invert(self):
        self.x = self.x
        self.y = self.y * -1

    def double(self):
        u_sq = -2
        m_numerator = 3 * self.x**2
        m_denominator = self.y * 2

        m_sq = m_numerator**2 / (m_denominator**2 * u_sq)
        new_x = m_sq - self.x * 2
        # Multiply m by u, which gets us to be able to substitute u^2 = -2 in the denominator.
        new_y = (m_numerator / m_denominator / u_sq) * (self.x * 3 - m_sq) - self.y

        self.x = new_x
        # y will be in terms of u.
        self.y = new_y

        return type(self)(self.x.value, self.y.value)


class G1:
    prime = 101

    def __init__(self, x, y):
        self.x = PrimeField(x)
        self.y = PrimeField(y)

    def invert(self):
        self.x = self.x
        self.y = self.y * -1

    def double(self):
        m = (3 * self.x**2) / (self.y * 2)
        new_x = m**2 - self.x * 2
        new_y = m * (self.x * 3 - m**2) - self.y

        self.x = new_x
        self.y = new_y

        return type(self)(self.x.value, self.y.value)


class Selector(Enum):
    L = 0
    R = 1
    O = 2
    M = 3
    C = 4


class Gate:
    """
    Represents a full PLONK gate:
    q_L * a + q_R + b + q_O * c + q_M * a * b + q_C = 0
    """

    def __init__(self, q_L, q_R, q_O, q_M, q_C):
        self.wires = {}
        self.wires[Selector.L] = q_L
        self.wires[Selector.R] = q_R
        self.wires[Selector.O] = q_O
        self.wires[Selector.M] = q_M
        self.wires[Selector.C] = q_C


def collect_selector_vector(gates, selector):
    vector = []
    for gate in gates:
        vector.append(gate.wires[selector])

    return vector


class TestPart1(unittest.TestCase):
    def test_collect_selector_vectors(self):
        gate_1 = Gate(0, 0, -1, 1, 0)
        gate_2 = Gate(0, 0, -1, 1, 0)
        gate_3 = Gate(0, 0, -1, 1, 0)
        gate_4 = Gate(1, 1, -1, 0, 0)

        q_L_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.L)
        q_R_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.R)
        q_O_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.O)
        q_M_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.M)
        q_C_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.C)

        expected_q_L = [0, 0, 0, 1]
        expected_q_R = [0, 0, 0, 1]
        expected_q_O = [-1, -1, -1, -1]
        expected_q_M = [1, 1, 1, 0]
        expected_q_C = [0, 0, 0, 0]

        self.assertEqual(q_L_vec, expected_q_L)
        self.assertEqual(q_R_vec, expected_q_R)
        self.assertEqual(q_O_vec, expected_q_O)
        self.assertEqual(q_M_vec, expected_q_M)
        self.assertEqual(q_C_vec, expected_q_C)

    def test_SRS(self):
        expected_srs = [
            G1(1, 2),
            G1(68, 74),
            G1(65, 98),
            G1(18, 49),
            G1(1, 99),
            G1(68, 27),
            G1(65, 3),
            G2(36, 31),
            G2(90, 82),
        ]

        self.assertEqual(StructuredReferenceString(), expected_srs)

    def test_doubling_and_invert(self):
        """
        Tests the doubling and inversion done in Part I for point doubling and inversion,
        up to 16G and -16G respectively.
        """
        p = G1(1, 2)
        inv_p = G1(1, 2)

        # -G = (68, 99)
        inv_p.invert()
        self.assertEqual(inv_p.x, 1)
        self.assertEqual(inv_p.y, 99)

        # 2G = (68, 74)
        p.double()
        self.assertEqual(p.x, 68)
        self.assertEqual(p.y, 74)

        # -2G = (68, 99)
        inv_p.double()
        self.assertEqual(inv_p.x, 68)
        self.assertEqual(inv_p.y, 27)

        # 4G = (65, 98)
        p.double()
        self.assertEqual(p.x, 65)
        self.assertEqual(p.y, 98)

        # -4G = (65, 3)
        inv_p.double()
        self.assertEqual(inv_p.x, 65)
        self.assertEqual(inv_p.y, 3)

        # 8G = (18, 49)
        p.double()
        self.assertEqual(p.x, 18)
        self.assertEqual(p.y, 49)

        # -8G = (18, 52)
        inv_p.double()
        self.assertEqual(inv_p.x, 18)
        self.assertEqual(inv_p.y, 52)

        # 16G = (1, 99)
        p.double()
        self.assertEqual(p.x, 1)
        self.assertEqual(p.y, 99)

        # -16G = (1, 2)
        inv_p.double()
        self.assertEqual(inv_p.x, 1)
        self.assertEqual(inv_p.y, 2)


if __name__ == "__main__":
    unittest.main()
