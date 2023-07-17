from enum import Enum
from field import Field17
from polynomial import Polynomial, interpolate_poly_degree_3
from curve import G1
from srs import StructuredReferenceString

# We use elliptic curve: y^2 = x^3 + 3

# Constants used to evaluate cosets.
# k1 is chosen such that it is not an element of H,
# k2 is chosen such that it is neither an element of H nor k1H.
k1 = 2
k2 = 3

n = 4


def preprocess_input():
    """Simulates input preprocessing, returning all selector polynomials."""
    roots_of_unity = expand_root_of_unity()
    gate_1 = Gate(0, 0, -1, 1, 0)  # a1 * b1 = c1
    gate_2 = Gate(0, 0, -1, 1, 0)  # a2 * b2 = c2
    gate_3 = Gate(0, 0, -1, 1, 0)  # a3 * b3 = c3
    gate_4 = Gate(1, 1, -1, 0, 0)  # a4 + b4 = c4

    q_L_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.L)
    q_R_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.R)
    q_O_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.O)
    q_M_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.M)
    q_C_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.C)

    q_L = interpolate_poly_degree_3(roots_of_unity, q_L_vec)
    q_R = interpolate_poly_degree_3(roots_of_unity, q_R_vec)
    q_O = interpolate_poly_degree_3(roots_of_unity, q_O_vec)
    q_M = interpolate_poly_degree_3(roots_of_unity, q_M_vec)
    q_C = interpolate_poly_degree_3(roots_of_unity, q_C_vec)

    return q_L, q_R, q_O, q_M, q_C


class Selector(Enum):
    """Selector wires for a PLONK gate."""
    L = 0
    R = 1
    O = 2
    M = 3
    C = 4


class Gate:
    """Returns a representation of a full PLONK gate.

    The representation is as such:
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
    """Collects the wires for a list of gates given a selector."""
    return [gate.wires[selector] for gate in gates]


def expand_root_of_unity(omega=4, modulus=17):
    """
    Given a root of unity (omega) and a modulo value of a chosen base field,
    calculate and return the roots of unity. Following the blog post,
    we use base field F_17 and choose omega = 4.
    """
    roots = [1]
    next_root = omega

    while next_root != 1:
        roots.append(next_root)
        next_root = roots[-1] * omega % modulus

    return roots


def evaluate_coset(roots_of_unity, constant, modulus=17):
    return [(root * constant) % modulus for root in roots_of_unity]


def make_copy_constraints():
    """Connect variables between gates by creating copy constraints (encoded as polynomials)."""
    roots_of_unity = expand_root_of_unity()
    H = roots_of_unity
    K1H = evaluate_coset(roots_of_unity, k1)
    K2H = evaluate_coset(roots_of_unity, k2)

    s_sigma1 = interpolate_poly_degree_3(
        roots_of_unity, [K1H[0], K1H[1], K1H[2], K2H[0]]
    )
    s_sigma2 = interpolate_poly_degree_3(roots_of_unity, [H[0], H[1], H[2], K2H[1]])
    s_sigma3 = interpolate_poly_degree_3(roots_of_unity, [H[3], K1H[3], K2H[3], K2H[2]])

    return s_sigma1, s_sigma2, s_sigma3


class Verifier:
    def __init__(self):
        """Inits a verifier with mock randomness simulated by dice rolls."""
        self.alpha = 15
        self.beta = 12
        self.gamma = 13
        self.zeta = 5

    def challenges(self):
        return self.alpha, self.beta, self.gamma, self.zeta

    def opening_challenge(self):
        """Returns the v found in Round 5."""
        return 12


class Prover:
    def __init__(self):
        pass

    def simulate_trusted_setup(self):
        """A mock trusted setup that returns a structured reference string.

        For our use case, this is a list of elliptic curve points
        calculated from our "secret" s.
        """
        self.s = 2
        self.srs = StructuredReferenceString(self.s)

    def eval_at_secret(self, vector):
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
        for coeff, point in zip(vector, self.srs):
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

    def encode_wire_polynomials(self):
        """Encode assignment vectors a, b and c for later use."""
        print("\n(Round 1) Encoding assignment vectors a, b, c for later use")
        # Generate random b_1, ... b6 from our base field (F17).
        # We pretend we rolled the die like in the blog to get:
        b = [7, 4, 11, 12, 16, 2]

        # Z_H is the polynomial that's zero on all elements of our
        # subgroup H. Since H is the 4th roots of unity in our field,
        # and x^4 = 1, x^4 - 1 = 0 is the zerofier.
        # We will express this in terms of coefficients, eg.
        # the value in i-th index is the coefficient of x^i.
        Z_H = Polynomial([-1, 0, 0, 0, 1])

        b1x_plus_b2 = Polynomial([4, 7])
        b3x_plus_b4 = Polynomial([12, 11])
        b5x_plus_b6 = Polynomial([2, 16])

        roots_of_unity = expand_root_of_unity()

        a = [3, 4, 5, 9]
        b = [3, 4, 5, 16]
        c = [9, 16, 25, 25]
        f_a = Polynomial(interpolate_poly_degree_3(roots_of_unity, a))
        f_b = Polynomial(interpolate_poly_degree_3(roots_of_unity, b))
        f_c = Polynomial(interpolate_poly_degree_3(roots_of_unity, c))

        # Compute wire polynomials
        a_x = Z_H * b1x_plus_b2 + f_a
        b_x = Z_H * b3x_plus_b4 + f_b
        c_x = Z_H * b5x_plus_b6 + f_c

        a_x_evaluated = self.eval_at_secret(a_x)
        b_x_evaluated = self.eval_at_secret(b_x)
        c_x_evaluated = self.eval_at_secret(c_x)

        print(f"a(x) = {a_x}, [a(s)] = {a_x_evaluated}")
        print(f"b(x) = {b_x}, [b(s)] = {b_x_evaluated}")
        print(f"c(x) = {c_x}, [c(s)] = {c_x_evaluated}\n")

        return (a_x, b_x, c_x, a_x_evaluated, b_x_evaluated, c_x_evaluated)

    def compute_permutation_poly(self, copy_constraints):
        """Commit to a single polynomial z that encodes all of the copy constraints."""
        print("(Round 2) Commit to z that encodes all copy constraints")
        s_sigma1, s_sigma2, s_sigma3 = copy_constraints

        # 1) Generate random b_1, ..., b_q in F_17
        # For simplicity, we use the same die results as in the blog.
        b_7 = 14
        b_8 = 11
        b_9 = 7

        # 2) Get challenges beta, gamma in F_17
        beta = 12
        gamma = 13

        # 3) Compute Z(x): Z(x) = (b_7 * x^2 + b_8 + b_9) * Z_H(x) + acc(x)
        # Z_H is the polynomial that's zero on all elements of our
        # subgroup H. Since H is the 4th roots of unity in our field,
        # and x^4 = 1, x^4 - 1 = 0 is the zerofier.
        # We will express this in terms of coefficients, eg.
        # the value in i-th index is the coefficient of x^i.
        Z_H = Polynomial([-1, 0, 0, 0, 1])

        a = [3, 4, 5, 9]
        b = [3, 4, 5, 16]
        c = [9, 16, 25, 25]

        # 4) output [Z(s)]
        omega = 4

        acc = [1]
        for j in range(1, 4):
            acc_value = 0
            for i in range(j):
                a_1 = Field17(a[i] + beta * (pow(omega, i)) + gamma)
                b_1 = Field17(b[i] + beta * k1 * pow(omega, i) + gamma)
                c_1 = Field17(c[i] + beta * k2 * pow(omega, i) + gamma)

                a_1_dem = Field17(a[i] + beta * s_sigma1.eval(pow(omega, i)) + gamma)
                b_1_dem = Field17(b[i] + beta * s_sigma2.eval(pow(omega, i)) + gamma)
                c_1_dem = Field17(c[i] + beta * s_sigma3.eval(pow(omega, i)) + gamma)

                num = a_1 * b_1 * c_1
                denom = a_1_dem * b_1_dem * c_1_dem

                if acc_value == 0:
                    acc_value = num / denom
                else:
                    acc_value *= num / denom

            acc.append(acc_value)

        roots_of_unity = expand_root_of_unity()
        acc_poly = interpolate_poly_degree_3(roots_of_unity, acc)

        z_poly = Polynomial([b_9, b_8, b_7]) * Z_H + acc_poly

        z_s_evaluated = self.eval_at_secret(Polynomial([z_poly.eval(self.s)]))

        print(f"z(x) = {z_poly}, [z(s)] = {z_s_evaluated}\n")
        return z_poly, z_s_evaluated

    def compute_quotient_poly(
        self,
        a_x,
        b_x,
        c_x,
        copy_constraints,
        challenges,
        preprocessed_input,
        pi_x,
        z_h,
        z_x,
    ):
        """Computes quotient polynomial t (degree = 3n + 5 for n gates).

        This polynomial encodes the majority of the info contained in our circuit
        and assignments all at once.
        """
        print("(Round 3) Computing quotient polynomial t(x) and its parts\n", z_h)

        ql_x, qr_x, qo_x, qm_x, qc_x = preprocessed_input
        alpha, beta, gamma, _ = challenges
        s_sigma1_x, s_sigma2_x, s_sigma3_x = copy_constraints

        t_1_x = a_x * b_x * qm_x + a_x * ql_x + b_x * qr_x + c_x * qo_x + pi_x + qc_x
        t_2_x = (
            Polynomial([alpha])
            * ((a_x + Polynomial([0, beta]) + gamma))
            * (b_x + (Polynomial([0, beta * k1]) + gamma))
            * (c_x + (Polynomial([0, beta * k2]) + gamma))
            * z_x
        )

        # This is a special case - we substitute x with omega * x,
        # i.e. at i-th position, it is omega^i.
        z_omega_x = Polynomial(
            [coefficient * pow(4, i) for i, coefficient in enumerate(z_x.coefficients)]
        )
        t_3_x = (
            Polynomial([alpha])
            * (a_x + s_sigma1_x * beta + gamma)
            * (b_x + s_sigma2_x * beta + gamma)
            * (c_x + s_sigma3_x * beta + gamma)
            * z_omega_x
        )

        roots_of_unity = expand_root_of_unity()
        L1_x = interpolate_poly_degree_3(roots_of_unity, [1, 0, 0, 0])
        t_4_x = (z_x - 1) * L1_x * Polynomial([pow(alpha, 2)])

        t_x = t_1_x + t_2_x - t_3_x + t_4_x

        t_x = t_x.synthetic_divide_similar(z_h)
        print("Quotient polynomials:\nt(x) =", t_x)

        # Break t(x) into three parts of degree n + 1 each
        t_lo = Polynomial(t_x[0:6])
        t_mid = Polynomial(t_x[6:12])
        t_hi = Polynomial(t_x[12:18])

        print(f"t_lo(x) = {t_lo}")
        print(f"t_mid(x) = {t_mid}")
        print(f"t_hi(x) = {t_hi}\n")

        return t_lo, t_mid, t_hi, t_x

    def compute_opening_evaluations(
        self,
        challenges,
        a_x,
        b_x,
        c_x,
        z_x,
    ):
        """Compute opening evaluations of our polynomials at zeta."""
        print("(Round 4) Compute opening evaluations")
        (_, _, _, zeta) = challenges

        a_hat = a_x.eval(zeta)
        b_hat = b_x.eval(zeta)
        c_hat = c_x.eval(zeta)

        omega = 4
        s_sigma1_hat = Polynomial([7, 13, 10, 6]).eval(zeta)
        s_sigma2_hat = Polynomial([4, 0, 13, 1]).eval(zeta)
        z_omega_hat = z_x.eval(zeta * omega)

        print(
            f"Evaluations:\na = {a_hat}, b = {b_hat}, c = {c_hat}\ns_sigma1 = {s_sigma1_hat}, s_sigma2 = {s_sigma2_hat}\nz_omega_hat = {z_omega_hat}\n"
        )

        return a_hat, b_hat, c_hat, s_sigma1_hat, s_sigma2_hat, z_omega_hat

    def compute_linearization_polynomial(
        self, challenges, opening_evaluations, z_x, preprocessed_input, pi_x
    ):
        """Compute linearization polynomial r(x)."""
        (
            a_hat,
            b_hat,
            c_hat,
            s_sigma1_hat,
            s_sigma2_hat,
            z_omega_hat,
        ) = opening_evaluations
        print("(Round 5) Compute linearization polynomial r(x)")
        alpha, beta, gamma, zeta = challenges
        ql_x, qr_x, qo_x, qm_x, qc_x = preprocessed_input

        s_sigma3 = Polynomial([6, 7, 3, 14])
        roots_of_unity = expand_root_of_unity()
        L1_zeta = interpolate_poly_degree_3(roots_of_unity, [1, 0, 0, 0]).eval(zeta)

        r_x = (
            (
                a_hat * b_hat * qm_x
                + a_hat * ql_x
                + b_hat * qr_x
                + c_hat * qo_x
                + pi_x
                + qc_x
            )
            + (
                alpha
                * (
                    (a_hat + beta * zeta + gamma)
                    * (b_hat + beta * k1 * zeta + gamma)
                    * (c_hat + beta * k2 * zeta + gamma)
                    * z_x
                )
                - (
                    alpha
                    * (a_hat + beta * s_sigma1_hat + gamma)
                    * (b_hat + beta * s_sigma2_hat + gamma)
                    * beta
                    * z_omega_hat
                    * s_sigma3
                )
            )
            + ((pow(alpha, 2) * L1_zeta) * z_x)
        )

        r_hat = r_x.eval(zeta)

        print(
            f"Linearization polynomial and evaluation:\nr(x) = {r_x}, r_hat = {r_hat}\n"
        )
        return r_x, r_hat

    def compute_opening_proof(
        self,
        zeta,
        v,
        copy_constraints,
        opening_evaluations,
        a_x,
        b_x,
        c_x,
        t_x,
        t_lo_x,
        t_mid_x,
        t_hi_x,
        z_x,
        r_x,
        r_hat,
    ):
        """Compute opening proofs."""
        (
            a_hat,
            b_hat,
            c_hat,
            s_sigma1_hat,
            s_sigma2_hat,
            z_omega_hat,
        ) = opening_evaluations
        (s_sigma1_x, s_sigma2_x, _) = copy_constraints

        w_zeta_x = (
            t_lo_x
            + pow(zeta, n + 2) * t_mid_x
            + pow(zeta, 2 * n + 4) * t_hi_x
            - t_x.eval(zeta)
            + v * (r_x - r_hat)
            + pow(v, 2) * (a_x - a_hat)
            + pow(v, 3) * (b_x - b_hat)
            + pow(v, 4) * (c_x - c_hat)
            + pow(v, 5) * (s_sigma1_x - s_sigma1_hat)
            + pow(v, 6) * (s_sigma2_x - s_sigma2_hat)
        )

        w_zeta_x = w_zeta_x.synthetic_divide(Polynomial([zeta, 1]))
        w_zeta_omega_x = z_x - z_omega_hat
        w_zeta_omega_x = w_zeta_omega_x.synthetic_divide(Polynomial([3, 1]))

        w_zeta_x_evaluated = self.eval_at_secret(w_zeta_x)
        w_zeta_omega_x_evaluated = self.eval_at_secret(w_zeta_omega_x)
        print(f"[Wz(x)] = {w_zeta_x_evaluated}")
        print(f"[Wzomega(x)] = {w_zeta_omega_x_evaluated}")

        return w_zeta_x_evaluated, w_zeta_omega_x_evaluated
