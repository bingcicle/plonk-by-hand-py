def Field(modulus):
    class FieldElement:
        def __init__(self, value):
            if isinstance(value, FieldElement):
                self.value = value.value % FieldElement.modulus
            else:
                self.value = value % FieldElement.modulus

        def __repr__(self):
            return f"{self.value}"

        def __add__(self, other):
            if isinstance(other, int):
                return FieldElement(self.value + other)

            return FieldElement(self.value + other.value)

        def __radd__(self, other):
            return self.__add__(other)

        def __sub__(self, other):
            if isinstance(other, int):
                return FieldElement(self.value - other)
            return FieldElement(self.value - other.value)

        def __rsub__(self, other):
            return self.__sub__(other)

        def __rmul__(self, other):
            return self.__mul__(other)

        def __mul__(self, other):
            from polynomial import Polynomial

            if isinstance(other, Polynomial):
                return other.__mul__(self)

            if isinstance(other, FieldElement):
                val = other.value
                return FieldElement(self.value * val)
            elif isinstance(other, int):
                val = other
                return FieldElement(self.value * val)

        def __pow__(self, other):
            if isinstance(other, FieldElement):
                value = other.value
            elif isinstance(other, int):
                value = other
            return FieldElement(self.value**value)

        def __truediv__(self, other):
            return self.__div__(other)

        def __rtruediv__(self, other):
            return self.__div__(other)

        def __neg__(self):
            return FieldElement(-self.value)

        def __div__(self, other):
            if isinstance(other, FieldElement):
                value = other
            elif isinstance(other, int):
                value = FieldElement(other)

            return FieldElement(self.value * value.inv())

        def __le__(self, other):
            if isinstance(other, int):
                return self.value <= other
            return self.value <= other.value

        def __lt__(self, other):
            if isinstance(other, int):
                return self.value < other
            return self.value < other.value

        def __gt__(self, other):
            if isinstance(other, int):
                return self.value > other
            return self.value > other.value

        def __ge__(self, other):
            if isinstance(other, int):
                return self.value >= other
            return self.value >= other.value

        def __eq__(self, other):
            if isinstance(other, FieldElement):
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

    FieldElement.modulus = modulus
    return FieldElement


Field17 = Field(17)
Field101 = Field(101)
