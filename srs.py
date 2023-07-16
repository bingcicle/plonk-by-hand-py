from curve import G1, G2


class StructuredReferenceString:
    """Represents a list of elliptic curve points parameterized
    by a randomly generated secret number s.
    """

    values = []

    def __init__(self, s):
        """Inits an SRS with 9 (n + 5) elements (since our chosen gate will
        have 4 (n) gates), using randomly chosen powers of s to construct
        this reference string.

        To follow the blog, we pretend to throw a die that results in s = 2,
        and calculate 6 powers of s to use as coefficients for our generators.
        """
        powers_of_s = [s * i for i in range(6)]

        G1(1, 2)
        first = [G1(1, 2)]
        for _ in range(len(powers_of_s)):
            first.append(first[-1].double())

        g2 = G2(36, 31)
        second = [G2(36, 31)]
        doubled = g2.double()
        second.append(doubled)

        self.values = first + second

    def __getitem__(self, i):
        return self.values[i]

    def __iter__(self):
        return iter(self.values)

    def __eq__(self, other):
        return all(
            a.x.value == b.x.value and a.y.value == b.y.value
            for (a, b) in zip(self, other)
        )

    def __repr__(self):
        return f"[{', '.join('(' + str(v.x.value) + ', ' + str(v.y.value) + ')' for v in self.values)}]"
