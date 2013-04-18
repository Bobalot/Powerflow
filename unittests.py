#!/usr/bin/python
#	Robert Williamson
#	2013

from powerflow import Netlist

import numpy as np
import time

import unittest
import nose


class TestFile1(unittest.TestCase):
	def runTest(self):
		# Already known solution
		node_list = ['n1', 'n2', 'n3', 'n4', 'n5']
		V = np.complex64(( 1.04496808+0.00816805j, 0.94438640-0.08014649j, 1.00000000+0.j, 0.94980016-0.11007433j, 0.96739431-0.07275j))
		phase = np.float64(( 0.00771616, -0.08476695, 0.0, -0.11545534, -0.07506073))

		P = np.float64((0.95, -0.8, 1.11453631, -0.7, -0.5))
		Q = np.float64((0.92777154, -0.65, -0.02975013, 0.1, -0.075))

		netlist = Netlist()
		# netlist.print_results()

		# Check values are close to what we expect
		
		try:
			netlist.load_from_file("testing/inputfile1.txt")
			netlist.run()

			assert (np.allclose(netlist.result.V, V, rtol=1e-3))
			assert (np.allclose(netlist.result.phase, phase, rtol=1e-3))
			assert (np.allclose(netlist.result.P, P, rtol=1e-3))
			assert (np.allclose(netlist.result.Q, Q, rtol=1e-3))
		finally:
			del netlist


class TestFile2(unittest.TestCase):
	def runTest(self):
		node_list = ['n1', 'n2', 'n3']
		V = np.complex64(( 1.02500000+0.0j, 1.00057074-0.03669173j, 1.02970624+0.0245978j))
		phase = np.complex64((0.0, -0.03665437, 0.02388363))

		P = np.complex64((1.00010597, -4.0, 3.0))
		Q = np.complex64((0.90512173, -2.0, 1.36935805))

		netlist = Netlist()

		
		# Check values are close to what we expect
		

		try:
			netlist.load_from_file("testing/inputfile2.txt")
			netlist.run()

			# print abs(netlist.result.V)- abs(V)

			print netlist.result.V
			print netlist.result.P


			# assert (np.allclose(netlist.result.V, V, rtol=1e-3, atol=1e-3))
			# assert (np.allclose(netlist.result.phase, phase, rtol=1e-3, atol=1e-3))
			# assert (np.allclose(netlist.result.P, P, rtol=1e-3, atol=1e-3))
			# assert (np.allclose(netlist.result.Q, Q, rtol=1e-3, atol=1e-3))
		finally:
			del netlist


if __name__=='__main__':
	
	# suite = unittest.TestLoader().loadTestsFromTestCase(TestFile1)
	# unittest.TextTestRunner(verbosity=2).run(suite)

	suite = unittest.TestLoader().loadTestsFromTestCase(TestFile2)
	unittest.TextTestRunner(verbosity=2).run(suite)



