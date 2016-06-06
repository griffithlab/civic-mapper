"""civic_mapper.py - Maps and reports variants found in CIViC using the CIViC API"""

import argparse
from flask import Flask, render_template
import json
import requests
requests.packages.urllib3.disable_warnings()
import re
import gzip


# Initializing parser to argparse.ArgumentParser
parser = argparse.ArgumentParser( description = __doc__)
# Specificing that the program requires the '--transvar' informer along with a transvar file
parser.add_argument('-transvar', '--transvar', required = False, help = 'Reports variants found both in the transvar file and CIViC')
# Initializing args to parse through all arguments
args = parser.parse_args()

# Saves the inputed transvar file as a readable variables
if args.transvar != None:
	transvar_output = open(args.transvar, 'r')


# This code accesses CIViC and creates dictionaries of variants with and without coordinates
# Variants with coordinates are organized as such; dictionary = {chromosome number: [list of variants in that chromosome]}
# Variants witout coordinates are organized as such; dictionary = {gene name: [list of variants in that gene]}
##################################################
def coordinate_sorter(start, stop, ascending_list):
	"""Maintains an ascending list by properly placing the next object in order when given the start and stop coordinates of the current object and list desired to be maintained.
		Note: The function works best starting with an empty list, or a list with a single entry, or an already organized ascending list

		Args:
			start = First coordinate
			stop = Last coordinate
			ascending_list = list of data already sorted by stop coordinates (least to greatest)

		Returns:
			The index location of where the current object should be inserted within the list
	"""
	if len(ascending_list) == 0:
		#print('1st IF')
		return 0
	else:
		for current_index in range(0, len(ascending_list)):

			if start == stop and start in range(int(ascending_list[current_index]['coordinates']['start']), int(ascending_list[current_index]['coordinates']['stop'])+1):
				#print('2nd IF')
				return current_index

			elif start == stop and stop >= int(ascending_list[current_index]['coordinates']['stop']) and current_index+1 == len(ascending_list):
				#print('3rd IF')
				return current_index+1

			elif stop >= int(ascending_list[current_index]['coordinates']['stop']) and current_index+1 == len(ascending_list):
				#print('4th IF')
				return current_index+1

# The civic_variants_with_coordinates_dictionary will hold all variants based on chromosome location (dictionary key) and exist within a list
civic_variants_with_coordinates_dictionary = {'X':[], 'Y':[]}
for chromosome_number in range(1,23):
	civic_variants_with_coordinates_dictionary[str(chromosome_number)] = []

raw_civic_variants_with_coordinates_list = []
# The civic_variants_without_coordinates_dictionary will hold all variants based on chromosome location (dictionary key) and exist within a list
civic_variants_without_coordinates_dictionary = {}
# Directly requesting a list of all variants from the CIViC API in JSON format
variants = requests.get('https://civic.genome.wustl.edu/api/variants?count=10000').json()['records']

# Iterating through the entire list of variants and parsing out variants without coordinates, while keeping count of all variants
for current_variant in range(0, len(variants)):
	#####Currently filtering out variants without representative trascripts because they most likely won't have coordinates
	"""NOTE: Variants that are [FUSIONS, EXPRESSION, UNDEREXPRESSION, OVEREXPRESSION, AMPLIFICATION, LOSS] are also being removed currently to be dealt with later COME BACK TO HANDLE THESE VARIANTS FOR MAPPING AND REPORTATION."""
	
	if variants[current_variant]['coordinates']['representative_transcript'] == None or variants[current_variant]['coordinates']['chromosome2'] != None or ('EXPRESSION' in variants[current_variant]['name']) or ('AMPLIFICATION' in variants[current_variant]['name']) or ('LOSS' in variants[current_variant]['name']) or ('3\' FUSION' in variants[current_variant]['name']) or ('LATION' in variants[current_variant]['name']):
		gene_name = variants[current_variant]['entrez_name']
		if gene_name in civic_variants_without_coordinates_dictionary:
			civic_variants_without_coordinates_dictionary[gene_name].append(variants[current_variant])
		else:
			civic_variants_without_coordinates_dictionary[gene_name] = [variants[current_variant]]
	else:
		raw_civic_variants_with_coordinates_list.append(variants[current_variant])

# Sorting all variants by stop coordinate in ascending order (least to greatest) 
raw_civic_variants_with_coordinates_list = sorted(raw_civic_variants_with_coordinates_list, key = lambda k: int(k['coordinates']['stop']))
for var in raw_civic_variants_with_coordinates_list:
	start = int(var['coordinates']['start'])
	stop = int(var['coordinates']['stop'])
	chr_key = var['coordinates']['chromosome']
	# Sorting nested and special cases of coordinates
	index_location = coordinate_sorter(start, stop, civic_variants_with_coordinates_dictionary[chr_key])
	# Inserting variants into civic_variants_with_coordinates_dictionary based on the index_location to maintain a sorted list
	civic_variants_with_coordinates_dictionary[chr_key].insert(index_location, var)

"""
for chromo in civic_variants_with_coordinates_dictionary:
	print('\n' + str(chromo))

	for var in civic_variants_with_coordinates_dictionary[chromo]:
		if int(var['coordinates']['start']) == int(var['coordinates']['stop']):
			print( var['entrez_name'] + ':p.' + var['name'] + '\t' + var['coordinates']['start'])
		else:
			print(var['entrez_name'] + ':p.' + var['name'] + '\t' + var['coordinates']['start'] + '\t' + var['coordinates']['stop'])

		#print(variants[current_variant]['entrez_name'] + ':p.' + variants[current_variant]['name'])
##################################################
#Extracts needed information from a transvar output file to map variants to CIViC and report findings
##################################################
"""
for line in transvar_output:

	line_list = line.strip().split()
	
	if line_list[0] != 'input':

		var_input = line_list[0]

		var_gene_name = line_list[3]

		var_var_name = ''
		temp_name = line_list[5].split('/')[-1]
		if 'p' in temp_name:
			var_var_name = line_list[5].split('/')[-1].split('.')[1]
			var_var_short_name = var_var_name[:-1]
		elif 'fs' in temp_name:
			var_var_name = line_list[5].split('/')[-1].split('.')[1]
			var_var_short_name = var_var_name.split('fs')[0][:-1]
		elif '.' == temp_name:
			var_var_name = None
		
		var_chr = ''
		for index in line_list[0].split(':')[0]:
			if index.isdigit() == True:
				var_chr += index

		variant_start = ''
		variant_stop = ''
		if '_' in line_list[0].split('.')[1]:
			variant_start = line_list[0].split('.')[1].split('_')[0]
			for index in line_list[0].split('_')[1]:
				if index.isdigit() == True:
					variant_stop += index
		else:
			for index in line_list[0].split('.')[1]:
				if index.isdigit() == True:
					variant_start += index
					variant_stop = variant_start

		var_ref_base = ''
		if '>' in line_list[0]:
			for index in line_list[0].split('.')[1].split('>')[0]:
				if index.isalpha() == True:
					var_ref_base = index
		elif 'ins' in line_list[0]:
			var_ref_base = None
		else:
			var_ref_base = line_list[7].split(';')[1].split('del')[1]

		var_var_base = ''
		if '>' in line_list[0]:
			var_var_base = line_list[0].split('>')[1]
		elif 'ins' in line_list[0]:
			var_var_base = line_list[0].split('ins')[1]
		else:
			var_var_base = None

		var_rep_trans = line_list[1]

		#########################input_genome_build = 

		passed_exact_match = False
		passed_soft_match = False
		passed_nested_match = False

		print_header_statement = False
		print_input_statement = False

		exact_match = ''

		soft_match = ''
		soft_matchs = []

		nested_match = ''
		nested_matchs = []

		for current_variant_in_dict in civic_variants_with_coordinates_dictionary[var_chr]:
			
			civic_start = int(current_variant_in_dict['coordinates']['start'])
			civic_stop = int(current_variant_in_dict['coordinates']['stop'])
			civic_ref_base = current_variant_in_dict['coordinates']['reference_bases']
			civic_var_base = current_variant_in_dict['coordinates']['variant_bases']
			civic_gene_name = current_variant_in_dict['entrez_name']
			civic_var_name = current_variant_in_dict['name']
			civic_rep_trans = current_variant_in_dict['coordinates']['representative_transcript'].split('.')[0]
			
			if var_rep_trans == civic_rep_trans and var_gene_name == civic_gene_name:

				## EXACT MATCH ##
				if int(variant_start) == civic_start and int(variant_stop) == civic_stop and var_ref_base == civic_ref_base and var_var_base == civic_var_base:
					passed_exact_match = True
					print_input_statement = True
					exact_match = current_variant_in_dict['entrez_name'] + ':' + current_variant_in_dict['name']

				## SOFT MATCH ##
				elif (passed_exact_match == False and passed_soft_match == False and var_var_short_name == civic_var_name and var_ref_base != civic_ref_base) or (passed_exact_match == False and passed_soft_match == False and var_ref_base != civic_ref_base and int(variant_start) in range(civic_start, civic_stop+1) and int(variant_stop) in range(civic_start, civic_stop+1)) :
					passed_soft_match = True
					print_input_statement = True
					soft_match = current_variant_in_dict['entrez_name'] + ':' + current_variant_in_dict['name']
					soft_matchs.append(soft_match)

				elif (passed_exact_match == True or passed_soft_match == True) and var_var_name != civic_var_name and var_ref_base != civic_ref_base and int(variant_start) in range(civic_start, civic_stop+1) and int(variant_stop) in range( civic_start, civic_stop+1):
					passed_nested_match = True
					print_header_statement = True
					nested_match = current_variant_in_dict['entrez_name'] + ':' + current_variant_in_dict['name']
					nested_matchs.append(nested_match)


		if print_input_statement == True:
			if var_var_name == None:
				print('\nVariant input information: ' + var_input + '\tAnnotated as: ' + var_gene_name + ':None')
			else:
				print('\nVariant input information: ' + var_input + '\tAnnotated as: ' + var_gene_name + ':' + var_var_name)
			
			if passed_exact_match == True:
				print('Exact match found to variants: ' + exact_match)
			if passed_soft_match == True:
				for curr_match in soft_matchs:
					print('Soft match to: ' + curr_match)
			if passed_nested_match == True:
				print('Variants below encompass the above variant found and might be of interest')
				for curr_match in nested_matchs:
					print('\t' + curr_match)

##################################################




