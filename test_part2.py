import unittest
from plonk import (
    Gate,
    expand_root_of_unity,
    Prover,
    Verifier,
    preprocess_input,
    collect_selector_vector,
    Selector,
    make_copy_constraints,
    evaluate_coset,
    k1,
    k2,
)
from curve import G1
from polynomial import Polynomial, interpolate_poly_degree_3


class TestPart2(unittest.TestCase):
    def test_interpolate(self):
        roots_of_unity = expand_root_of_unity()
        B = [3, 4, 5, 9]  # Our expected output wires
        f_a = interpolate_poly_degree_3(roots_of_unity, B)
        expected_f_a = Polynomial([1, 13, 3, 3])

        self.assertEqual(f_a, expected_f_a)

    def test_interpolate_circuit_vectors(self):
        gate_1 = Gate(0, 0, -1, 1, 0)
        gate_2 = Gate(0, 0, -1, 1, 0)
        gate_3 = Gate(0, 0, -1, 1, 0)
        gate_4 = Gate(1, 1, -1, 0, 0)

        q_L_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.L)
        q_R_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.R)
        q_O_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.O)
        q_M_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.M)
        q_C_vec = collect_selector_vector([gate_1, gate_2, gate_3, gate_4], Selector.C)

        roots_of_unity = expand_root_of_unity()

        expected_f_a = Polynomial([1, 13, 3, 3])
        expected_f_b = Polynomial([7, 3, 14, 13])
        expected_f_c = Polynomial([6, 5, 11, 4])
        expected_q_L = Polynomial([13, 1, 4, 16])
        expected_q_R = Polynomial([13, 1, 4, 16])
        expected_q_O = Polynomial([16, 0, 0, 0])
        expected_q_M = Polynomial([5, 16, 13, 1])
        expected_q_C = Polynomial([0, 0, 0, 0])

        a = Polynomial([3, 4, 5, 9])
        b = Polynomial([3, 4, 5, 16])
        c = Polynomial([9, 16, 25, 25])
        f_a = interpolate_poly_degree_3(roots_of_unity, a)
        f_b = interpolate_poly_degree_3(roots_of_unity, b)
        f_c = interpolate_poly_degree_3(roots_of_unity, c)
        q_L = interpolate_poly_degree_3(roots_of_unity, q_L_vec)
        q_R = interpolate_poly_degree_3(roots_of_unity, q_R_vec)
        q_O = interpolate_poly_degree_3(roots_of_unity, q_O_vec)
        q_M = interpolate_poly_degree_3(roots_of_unity, q_M_vec)
        q_C = interpolate_poly_degree_3(roots_of_unity, q_C_vec)

        self.assertTrue((f_a == expected_f_a))
        self.assertTrue((f_b == expected_f_b))
        self.assertTrue((f_c == expected_f_c))
        self.assertTrue((q_L == expected_q_L))
        self.assertTrue((q_R == expected_q_R))
        self.assertTrue((q_O == expected_q_O))
        self.assertTrue((q_M == expected_q_M))
        self.assertTrue((q_C == expected_q_C))

    def test_root_of_unity_with_cosets(self):
        H = [1, 4, 16, 13]
        K1H = [2, 8, 15, 9]
        K2H = [3, 12, 14, 5]

        roots_of_unity = expand_root_of_unity()
        self.assertEqual(roots_of_unity, H)
        self.assertEqual(evaluate_coset(roots_of_unity, k1), K1H)
        self.assertEqual(evaluate_coset(roots_of_unity, k2), K2H)

    def test_round_1(self):
        prover = Prover()
        prover.simulate_trusted_setup()
        (
            _,
            _,
            _,
            a_x_evaluated,
            b_x_evaluated,
            c_x_evaluated,
        ) = prover.encode_wire_polynomials()

        # Assert that after round 1, we get the correct points as shown
        # in the blog.
        self.assertEqual(a_x_evaluated, G1(91, 66))
        self.assertEqual(b_x_evaluated, G1(26, 45))
        self.assertEqual(c_x_evaluated, G1(91, 35))

    def test_round_2(self):
        prover = Prover()
        prover.simulate_trusted_setup()

        copy_constraints = make_copy_constraints()

        z_x, z_s_evaluated = prover.compute_permutation_poly(copy_constraints)

        self.assertEqual(z_s_evaluated, G1(32, 59))

    def test_round_3(self):
        prover = Prover()
        prover.simulate_trusted_setup()

        verifier = Verifier()
        challenges = verifier.challenges()
        preprocessed_input = preprocess_input()
        copy_constraints = make_copy_constraints()

        (
            a_x,
            b_x,
            c_x,
            _,
            _,
            _,
        ) = prover.encode_wire_polynomials()

        z_h = Polynomial([-1, 0, 0, 0, 1])

        z_x, _ = prover.compute_permutation_poly(copy_constraints)
        t_lo, t_mid, t_hi, t_x = prover.compute_quotient_poly(
            a_x,
            b_x,
            c_x,
            copy_constraints,
            challenges,
            preprocessed_input,
            Polynomial([0]),
            z_h,
            z_x,
        )
        t_lo_evaluated = prover.eval_at_secret(t_lo)
        t_mid_evaluated = prover.eval_at_secret(t_mid)
        t_hi_evaluated = prover.eval_at_secret(t_hi)

        self.assertEqual(t_lo_evaluated, G1(12, 32))
        self.assertEqual(t_mid_evaluated, G1(26, 45))
        self.assertEqual(t_hi_evaluated, G1(91, 66))

    def test_prover(self):
        prover = Prover()
        prover.simulate_trusted_setup()

        verifier = Verifier()

        preprocessed_input = preprocess_input()
        (
            a_x,
            b_x,
            c_x,
            a_x_evaluated,
            b_x_evaluated,
            c_x_evaluated,
        ) = prover.encode_wire_polynomials()
        copy_constraints = make_copy_constraints()

        z_h = Polynomial([-1, 0, 0, 0, 1])

        # chosen by verifier
        challenges = verifier.challenges()
        (_, _, _, zeta) = challenges

        z_x, z_x_evaluated = prover.compute_permutation_poly(copy_constraints)
        t_lo_x, t_mid_x, t_hi_x, t_x = prover.compute_quotient_poly(
            a_x,
            b_x,
            c_x,
            copy_constraints,
            challenges,
            preprocessed_input,
            Polynomial([0]),
            z_h,
            z_x,
        )

        t_lo_evaluated = prover.eval_at_secret(t_lo_x)
        t_mid_evaluated = prover.eval_at_secret(t_mid_x)
        t_hi_evaluated = prover.eval_at_secret(t_hi_x)

        pi_x = Polynomial([0])

        # Round 4
        opening_evaluations = prover.compute_opening_evaluations(
            challenges,
            a_x,
            b_x,
            c_x,
            z_x,
        )

        r_x, r_hat = prover.compute_linearization_polynomial(
            challenges, opening_evaluations, z_x, preprocessed_input, pi_x
        )
        v = verifier.opening_challenge()
        (
            a_hat,
            b_hat,
            c_hat,
            s_sigma1_hat,
            s_sigma2_hat,
            z_omega_hat,
        ) = opening_evaluations
        w_zeta_x_evaluated, w_zeta_omega_x_evaluated = prover.compute_opening_proof(
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
        )

        expected_proof = (
            G1(91, 66),
            G1(26, 45),
            G1(91, 35),
            G1(32, 59),
            G1(12, 32),
            G1(26, 45),
            G1(91, 66),
            G1(91, 35),
            G1(65, 98),
            15,
            13,
            5,
            1,
            12,
            15,
            15,
        )
        proof = (
            a_x_evaluated,
            b_x_evaluated,
            c_x_evaluated,
            z_x_evaluated,
            t_lo_evaluated,
            t_mid_evaluated,
            t_hi_evaluated,
            w_zeta_x_evaluated,
            w_zeta_omega_x_evaluated,
            a_hat,
            b_hat,
            c_hat,
            s_sigma1_hat,
            s_sigma2_hat,
            r_hat,
            z_omega_hat,
        )

        print(f"proof = {proof}")
        self.assertEqual(proof, expected_proof)
        print(f"proof generated correctly!")
