import unittest
from plonk import Gate, collect_selector_vector, StructuredReferenceString, Selector
from curve import G1, G2, generate_g1_subgroup


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

        self.assertEqual(StructuredReferenceString(s=2), expected_srs)

    def test_doubling_and_invert(self):
        """
        Tests the doubling and inversion done in Part 1 for point doubling and inversion,
        up to 16G and -16G respectively.
        """
        p = G1(1, 2)
        inv_p = G1(1, 2)

        # -G = (68, 99)
        inv_p.invert()
        self.assertEqual(inv_p.x, 1)
        self.assertEqual(inv_p.y, 99)

        # 2G = (68, 74)
        p = p.double()
        self.assertEqual(p.x, 68)
        self.assertEqual(p.y, 74)

        # -2G = (68, 99)
        inv_p = inv_p.double()
        self.assertEqual(inv_p.x, 68)
        self.assertEqual(inv_p.y, 27)

        # 4G = (65, 98)
        p = p.double()
        self.assertEqual(p.x, 65)
        self.assertEqual(p.y, 98)

        # -4G = (65, 3)
        inv_p = inv_p.double()
        self.assertEqual(inv_p.x, 65)
        self.assertEqual(inv_p.y, 3)

        # 8G = (18, 49)
        p = p.double()
        self.assertEqual(p.x, 18)
        self.assertEqual(p.y, 49)

        # -8G = (18, 52)
        inv_p = inv_p.double()
        self.assertEqual(inv_p.x, 18)
        self.assertEqual(inv_p.y, 52)

        # 16G = (1, 99)
        p = p.double()
        self.assertEqual(p.x, 1)
        self.assertEqual(p.y, 99)

        # -16G = (1, 2)
        inv_p = inv_p.double()
        self.assertEqual(inv_p.x, 1)
        self.assertEqual(inv_p.y, 2)

    def test_generate_g1_subgroup(self):
        subgroup = generate_g1_subgroup()
        expected_subgroup = [
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
        self.assertEqual(subgroup, expected_subgroup)
