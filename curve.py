from field import Field101


def generate_g1_subgroup():
    p = G1(1, 2)  # G1
    subgroup = [G1(1, 2), p.double()]

    # base field 17 - 1 - 2 elements = 14
    for _ in range(14):
        subgroup.append(subgroup[-1] + p)

    return subgroup


class G1:
    def __init__(self, x, y):
        self.x = Field101(x)
        self.y = Field101(y)

    def __repr__(self):
        return f"[{self.x.value}, {self.y.value}]"

    def __eq__(self, other):
        return self.x.value == other.x.value and self.y.value == other.y.value

    def invert(self):
        self.x = self.x
        self.y = self.y * -1

    def __add__(self, other):
        """Elliptic curve addition."""
        lambda_ = (other.y - self.y) / (other.x - self.x)
        x_r = pow(lambda_, 2) - self.x - other.x
        y_r = lambda_ * (self.x - x_r) - self.y
        return G1(x_r, y_r)

    def double(self):
        """Elliptic curve point doubling."""
        m = (3 * pow(self.x, 2)) / (self.y * 2)
        new_x = pow(m, 2) - self.x * 2
        new_y = m * (self.x * 3 - pow(m, 2)) - self.y

        return type(self)(new_x.value, new_y.value)


class G2:
    def __init__(self, x, y):
        self.x = Field101(x)
        self.y = Field101(y)

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
