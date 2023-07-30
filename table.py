from curve import G1
from field import Field17

class Table:
    """Represents the table containing all points as multiples of our generator G1.

    This is used to 'cheat' in our computation by looking up the results of our
    computations based on the multiples. Note that this cheat probably only makes sense
    when done by hand - we can definitely calculate it properly here with Python,
    but we stay true to the blog post here.
    """
    cheat_table = [
        G1(1, 2),
        G1(68, 74),
        G1(26, 45),
        G1(65, 98),
        G1(12, 32),
        G1(32, 42),
        G1(91, 35),
        G1(18, 49),
        G1(18, 52),
        G1(91, 66),
        G1(32, 59),
        G1(12, 69),
        G1(65, 3),
        G1(26, 56),
        G1(68, 27),
        G1(1, 99),
    ]

    def __init__(self):
        pass

    def eval_at_secret(self, vector, srs):
        """
        Evaluate elliptic curve points which represent
        polynomials evaluated at the "secret" s.
        """
        cheat_table = [
            G1(1, 2),
            G1(68, 74),
            G1(26, 45),
            G1(65, 98),
            G1(12, 32),
            G1(32, 42),
            G1(91, 35),
            G1(18, 49),
            G1(18, 52),
            G1(91, 66),
            G1(32, 59),
            G1(12, 69),
            G1(65, 3),
            G1(26, 56),
            G1(68, 27),
            G1(1, 99),
        ]

        # For each point, what we want to do is split up the doubles, i.e.
        # 13G = 8G + 4G + 1G.
        # Then based on the index i, power it by the same number as the coefficients, mod the prime.
        # eg. for 8G, if G is (68, 74), it's in the 2nd spot so double of 2 = 4.
        points = []
        for coeff, point in zip(vector, srs):
            # +1 to work in terms of mod 17, our chosen field.
            idx = cheat_table.index(point) + 1
            additions = []
            doubles = []

            # For each term, count the coefficient for the number of point doublings.
            while coeff > 1:
                num_doubles = 0
                powers_of_two = pow(2, num_doubles + 1)
                while powers_of_two <= coeff:
                    num_doubles += 1
                    powers_of_two = pow(2, num_doubles + 1)
                doubles.append(num_doubles)
                coeff -= pow(2, num_doubles)

            # If the coefficient is odd, we will be left with term, so we add it to the
            # final list to be summed.
            if coeff == 1:
                points.append(idx)

            # Once we've expressed this term in terms of doubles, double all together.
            for num in doubles:
                doubled = idx
                while num != 0:
                    doubled = doubled * 2
                    num -= 1

                additions.append(doubled)

            while len(additions) != 0:
                a = additions.pop()
                if len(additions) > 0:
                    b = additions.pop()
                    idx = a + b
                    additions.append(idx)
                else:
                    points.append(Field17(a))

        # All points now have coefficient = 1. Sum all points.
        final_point = points[0]
        for point in points[1 : len(points)]:
            final_point = final_point + point

        # Before we return, we need to subtract by 1 to index the correct term in the 0-15 table.
        return cheat_table[(final_point - 1).value]
