#!/usr/bin/python
#	Robert Williamson
#	2013

import numpy as np
import math
import string
import time
import cmath
import argparse

from collections import namedtuple
from tools import rect, cphase


class Element:
	# Reference numbers, so we can see the type of a bus and act on it
	SLACKTYPE=1
	GENTYPE=2
	LOADTYPE=3
	TRANSFORMERTYPE=4
	LINETYPE=5

	def __init__(self, name, ntype, node1, node2=0, **kwargs):
		self.name = name
		self.ntype = ntype

		self.node1 = node1
		self.node2 = node2

		self.calculated = dict()

		# Add node values
		self.V = np.float64( kwargs.get('V', 0) )
		self.phase = np.float64( kwargs.get('phase', 0) )
		self.P = np.float64( kwargs.get('P', 0) )
		self.Q = np.float64( kwargs.get('Q', 0) )

		# Add line values
		self.real = np.float64( kwargs.get('real', 0) )
		self.reactive = np.float64( kwargs.get('reactive', 0) )

Result = namedtuple('Result', ['node_list', 'node_types', 'I', 'Y', 'V', 'P', 'Q', 'phase', 'full_list', 'slacknodename' ], verbose=False)

class Netlist:

	device_list = set()
	device_list_names = set()

	# node_list = set()
	node_list_names = set()

	line_list = set()
	line_list_names = set()

	def __init__(self, verbose=False, slow_verbose=False):
		self.verbose = verbose
		self.slow_verbose = slow_verbose

	def count_node_types(self):
		slack_count = sum(1 for x in self.device_list if x.type == SLACKTYPE)
		gen_count = sum(1 for x in self.device_list if x.type == GENTYPE)
		load_count = sum(1 for x in self.device_list if x.type == LOADTYPE)

		return (slack_count, gen_count, load_count)

	def count_lines(self):
		return sum(1 for x in self.device_list if x.type == LINETYPE)

	def add_device(self, device_name ,device_type, node1, node2=0, **kwargs):

		if device_name not in self.device_list_names:
			self.device_list_names.add(device_name)

			new_device = Element(device_name, device_type, node1, node2, **dict(kwargs) )
			self.device_list.add(new_device)

			# Try and add the referenced nodes
			self.node_list_names.add(node1)
			self.node_list_names.add(node2)

	def add_line(self, line_name, node1, node2, **kwargs):
		kwreal = np.float64( kwargs.get('real', 0.0))
		kwreactive = np.float64( kwargs.get('reactive', 0.0))

		if line_name not in self.line_list_names:
			self.line_list_names.add(line_name)

			new_line = Element(line_name, Element.LINETYPE, node1, node2, real=kwreal, reactive=kwreactive )
			new_line.impedance = np.complex64((kwreal + (1j * kwreactive) ))
						
			self.line_list.add(new_line)

			# Try and add the referenced nodes
			self.node_list_names.add(node1)
			self.node_list_names.add(node2)


	def load_from_file(self, filename):
		# Check file exists
		try:
   			with open(filename) as f: pass
		except IOError as e:
   			print 'File does not exist'
   			exit()

		# Open file
		lines = [line.strip() for line in open(filename)]

		for line in lines:
			line2 = line.split()

			# Ignore blank lines and commented # lines
			if line2 != [] and line2[0].strip()[0]!="#":

				devicename = str(line2[0]).lower()
				node1 = str(line2[1]) 
				node2 = str(line2[2])

				# Get arguments printed in file
				arg1 = np.float64(str(line2[3]))
				arg2 = np.float64(str(line2[4]))


				if "z" in devicename:
					self.add_line(devicename, node1, node2, real=arg1, reactive=arg2 )

				elif "gen" in devicename:
					self.add_device(devicename, Element.GENTYPE, node1, node2, V=arg1, P=arg2)

				elif "slack" in devicename:
					self.add_device(devicename, Element.SLACKTYPE, node1, node2, V=arg1, phase=arg2)

				elif "load" in devicename:
					self.add_device(devicename, Element.LOADTYPE, node1, node2, P=arg1, Q=arg2)

				# devicetype = str(device[0]).lower()
				# if we have an impedance line
				# if "z" in devicetype[0]:
				# add nodes
				# self.node_list_names.add(line2[1])
				# self.node_list_names.add(line2[2])


	def netlist_sanity_check(self):
		nodelist = set(self.nodelist)
		is_sane = True

		# Check there is at least 1 slack node defined
		if self.count_lines[0]==0:
			print "No Slack nodes"
			is_sane = False

		# Check there are not multiple slack nodes
		if self.count_lines[0]>1:
			print "Multiple Slack nodes"
			is_sane = False

		return is_sane

	def run(self, rtol=1e-5, max_iterations=100000):
		# Copy the nodes to here so nothing gets overwriten
		device_list = list(self.device_list.copy())
		device_list_names = self.device_list_names.copy()

		node_list = list(self.node_list_names.copy())

		line_list = list(self.line_list.copy())
		line_list_names = self.line_list_names.copy()

		slacknodename = ''

		# Remove '0', the ground node from our node list
		if '0' in node_list: 
			node_list.remove("0")

		# count number of elements in the list
		numberofnodes = len(node_list)

		I = np.zeros((numberofnodes), dtype=np.complex64)

		Y = np.zeros((numberofnodes,numberofnodes), dtype=np.complex64)
		V = np.zeros((numberofnodes), dtype=np.complex64)

		P = np.zeros((numberofnodes), dtype=np.float64)
		Q = np.zeros((numberofnodes), dtype=np.float64)
		phase = np.zeros((numberofnodes), dtype=np.float64)

		node_types = np.zeros((numberofnodes),dtype=np.int32)

		# Set initial voltage on all nodes to 1+0j
		V[:] = 1+0j

		full_list = device_list
		full_list.extend(line_list)

		for node in node_list:
			for device in full_list:
				if device.node1 == node or device.node2 == node:
					devicetype = device.ntype

					# If we have an impedance line
					if device.ntype == Element.LINETYPE:
						impedance = device.impedance
						admittance = 1/impedance

						# If neither of the connected nodes are the ground node then build the admittance matrix normally
						if device.node1 !='0' and device.node2!='0':
							if node== device.node1:
								currentnode = node_list.index(device.node1)
								transfernode = node_list.index(device.node2)

							else:
								currentnode = node_list.index(device.node2)
								transfernode = node_list.index(device.node1)

							Y[currentnode][currentnode] += admittance
							Y[currentnode][transfernode] -= admittance

						# Check if either of the nodes are connected to ground 0
						elif device.node1 == '0':
							currentnode = node_list.index(device.node2)
							transfernode = currentnode
							Y[currentnode][transfernode] += admittance

						elif device.node2 == '0':
							currentnode = node_list.index(device.node1)
							transfernode = currentnode
							Y[currentnode][transfernode] += admittance

					
					elif device.ntype == Element.SLACKTYPE:
						idx = node_list.index(device.node1)

						V[idx] = rect( device.V, device.phase )
						phase[idx] = device.phase
						node_types[idx] = Element.SLACKTYPE 
						slacknodename = device.node1

					elif device.ntype == Element.GENTYPE:
						idx = node_list.index(device.node1)

						V[idx] = device.V
						P[idx] += device.P
						node_types[idx] = Element.GENTYPE

					elif device.ntype == Element.LOADTYPE:
						idx = node_list.index(device.node1)

						P[idx] -= device.P
						Q[idx] -= device.Q
						node_types[idx] = Element.LOADTYPE

		# Run calculations here
		m=0

		while m < max_iterations:


			Vold = V.copy()
			Qold = Q.copy()
			
			# start on the first node, node 0 is not the ground node, it is the first node in the list
			i=0
			while i<numberofnodes:

				# Check for generator node
				if node_types[i] == Element.GENTYPE:
					# sum1 is the sum of Y[i][j]*V[j], where i!=j. jc is used here instead of j, as j is used as the imaginary unit
					sum1 = np.complex64(0.0)
					for jc in range(0, numberofnodes):
						if jc!=i:
							sum1 += Y[jc,i]*V[jc]

					# Calculate new reactive power
					Q[i] = -(((V[i].conjugate()*( sum1 + Y[i][i]*V[i]) ).imag))

					# Calculate new phase
					S = np.complex64( P[i] - (1j * Q[i] ) )

					phase[i] = cphase( (1/Y[i][i]) * ((S / V[i].conjugate())-sum1))

					# Update V[i], so the magnitude is still the same and we have the new phase
					V[i] = rect( np.abs(V[i]), phase[i])


				# Check for load node
				elif node_types[i] == Element.LOADTYPE:

					# Sum of Y[i][j]*V[j], where i!=j
					sum1 = np.complex64(0.0)
					for jc in range(0, numberofnodes):
						if jc!=i:
							sum1 += Y[jc,i]*V[jc]

					S = np.complex64( P[i] - (1j*Q[i]) )
					V[i] = (1/Y[i][i]) * ((S / V[i].conjugate())-sum1)

					phase[i] = cphase(V[i])


				i+=1

			m+=1

			# Print iterations when in verbose mode
			if self.verbose:
				decimalplaces = 5

				print "Iteration: "+str(m)
				print "\n"+str(m)+" Node list"
				print node_list

				print "\n"+str(m)+" V"
				print np.around(V.transpose(), decimalplaces)

				print "\n"+str(m)+" abs(V)"
				print np.around(abs(V).transpose(), decimalplaces)

				print "\n"+str(m)+" phase (degrees)"
				print np.around(phase.transpose(), decimalplaces)*(180/np.pi)

				print "\n"+str(m)+" P"
				print np.around(P.transpose(), decimalplaces)

				print "\n"+str(m)+" Q"
				print np.around(Q.transpose(), decimalplaces)
			
				if self.slow_verbose:
					time.sleep(0.5)

			# If the deltaV/V for all its elements is less than the accuracy in the arguments,
			# usually 0.1%, then break the iteration loop

			# Ensure at least 1 iteration has occured
			if m>1:
				if all(abs(V[ V!=0 ]-Vold[ Vold!=0 ])/(Vold[ Vold!=0 ]) < rtol) and all(abs((Q[ Q!=0] - Qold[ Qold!= 0 ]))/Qold[ Qold!=0 ]  < rtol):
					# print abs((V-Vold)/Vold)
					if self.verbose:
						print "converged to within "+str(rtol*100)+"%  after "+str(m)+" iterations"
					break



		# Update the slack node 
		# Calculate P and Q on the slack bus
		Stotal = np.complex64(0+0j)

		# index of the slack node
		sidx = node_list.index(slacknodename)

		for device in device_list:
			if device.ntype == Element.LINETYPE:

				# See if one of these lines is connected to the slack node
				if device.node1==slacknodename or device.node2==slacknodename:
					
					# set snode to which node is connected and rnode to the remote node across the line
					if device.node1 == slacknodename:
						ridx = node_list.index(device.node2)

					elif device.node2 == slacknodename:
						ridx = node_list.index(device.node1)

					# Voltage at the slack node and the remote node
					vslacknode = V[sidx]
					vremote = V[ridx]

					admittance = 1/device.impedance

					lineflow = vslacknode*(vslacknode.conjugate() - vremote.conjugate())*admittance.conjugate()

					Stotal += lineflow

		# Update the P and Q on the slack node
		P[sidx] = Stotal.real
		Q[sidx] = Stotal.imag


		# Calculate the line flows



		result = Result(node_list=node_list,
				node_types=node_types,
				I=I,
				Y=Y,
				V=V,
				P=P,
				Q=Q,
				phase=phase,
				full_list=full_list,
				slacknodename=slacknodename	)

		self.result = result
		return result


	def get_line_flow(self, device_name):
		for device in self.result.full_list:
			if device.name == device_name and device.ntype==Element.LINETYPE:
				break

		V1 = complex(self.result.V[self.result.node_list.index(device.node1)])
		V2 = complex(self.result.V[self.result.node_list.index(device.node2)])
		
		P1 = complex(self.result.P[self.result.node_list.index(device.node1)])
		P2 = complex(self.result.P[self.result.node_list.index(device.node2)])
		
		Q1 = complex(self.result.Q[self.result.node_list.index(device.node1)])
		Q2 = complex(self.result.Q[self.result.node_list.index(device.node2)])
		
		admittance = 1 / device.impedance

		lineflow12 = V1*(V1.conjugate() - V2.conjugate())*admittance.conjugate()
		lineflow21 = V2*(V2.conjugate() - V1.conjugate())*admittance.conjugate()
		
		lineloss = lineflow12+lineflow21

		
		returndict = {'name':device.name,
					  'node1': device.node1,
					  'node2': device.node2,
					  'impedance': device.impedance,
					  'flow1to2': lineflow12,
					  'flow2to1': lineflow21,
					  'lineloss': lineloss}

		return returndict


	def print_results(self, decimalplaces=10):

		print "\nAdmittance Matrix"
		print np.around(self.result.Y, decimalplaces)

		print "\nNode list"
		print self.result.node_list
		# print nodetypes.transpose()

		print "\nV"
		print np.around(self.result.V, decimalplaces)

		print "\nabs(V)"
		print np.around(np.abs(self.result.V), decimalplaces)

		print "\nphase (degrees)"
		print np.around(self.result.phase*(180/math.pi), decimalplaces)

		print "\nP"
		print np.around(self.result.P, decimalplaces)

		print "\nQ"
		print np.around(self.result.Q, decimalplaces)

		print "\n"

		Plosses = np.sum(self.result.P)
		Qlosses = np.sum(self.result.Q)

		print "Total line losses (Real)"
		print np.around(Plosses, decimalplaces)
		print "Total line losses (Reactive)"
		print np.around(Qlosses, decimalplaces)
		print "\n"

		# Print the line flows and losses
		for device in self.result.full_list:
			if device.ntype == Element.LINETYPE:
				# V1 = complex(self.result.V[self.result.node_list.index(device.node1)])
				# V2 = complex(self.result.V[self.result.node_list.index(device.node2)])
				
				# P1 = complex(self.result.P[self.result.node_list.index(device.node1)])
				# P2 = complex(self.result.P[self.result.node_list.index(device.node2)])
				
				# Q1 = complex(self.result.Q[self.result.node_list.index(device.node1)])
				# Q2 = complex(self.result.Q[self.result.node_list.index(device.node2)])
				
				# admittance = 1 / device.impedance

				# lineflow12 = V1*(V1.conjugate() - V2.conjugate())*admittance.conjugate()
				# lineflow21 = V2*(V2.conjugate() - V1.conjugate())*admittance.conjugate()
				
				# lineloss = lineflow12+lineflow21

				line_flow = self.get_line_flow(device.name)


				print "From "+str(line_flow['node1'])+" to "+str(line_flow['node2'])+" through "+str(line_flow['name'])
				print "\tline flow leaving "+str(line_flow['node1'])+": "+str(np.around(line_flow['flow1to2'], decimalplaces))
				print "\tline flow leaving "+str(line_flow['node2'])+": "+str(np.around(line_flow['flow2to1'], decimalplaces))
				
				print "\tline loss: "+str(np.around(line_flow['lineloss'], decimalplaces))+"\n"









# Run the main() function
if __name__ == "__main__":
    # main()

	netlist = Netlist()
	netlist.load_from_file("testing/inputfile1.txt")
	netlist.run()
	netlist.print_results()


	netlist = Netlist()
	netlist.load_from_file("testing/inputfile2.txt")	
	netlist.run()
	netlist.print_results()

