"""civic_mapper.py - Maps and reports variants found in CIViC using the CIViC API"""

import argparse, json, requests, subprocess, transvar
from flask import Flask, render_template
requests.packages.urllib3.disable_warnings()



# Initializing parser to argparse.ArgumentParser
parser = argparse.ArgumentParser( description = __doc__)
# Specificing that the program requires the '--vcf' informer along with a vcf file
parser.add_argument('-vcf', '--vcf', required = False, help = 'Reports variants found both in the VCF file and CIViC')
# Initializing args to parse through all arguments
args = parser.parse_args()

# Saves the inputed VCF file as a readable variables
if args.vcf != None:
	#vcf_output = open(args.vcf, 'r')
	transvar_output = subprocess.check_output("transvar ganno --vcf " + args.vcf + " --ensembl --seqmax -1", shell = True)



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


# This code accesses CIViC and creates dictionaries of variants with and without coordinates
# Variants with coordinates are organized as such; dictionary = {chromosome number: [list of variants in that chromosome]}
# Variants witout coordinates are organized as such; dictionary = {gene name: [list of variants in that gene]}
##################################################
# The civic_variants_with_coordinates_dictionary will hold all variants based on chromosome location (dictionary key) and exist within a list
civic_variants_with_coordinates_dictionary = {'X':[], 'Y':[]}
for chromosome_number in range(1,23):
	civic_variants_with_coordinates_dictionary[str(chromosome_number)] = []

raw_civic_variants_with_coordinates_list = []
civic_variants_without_coordinates_dictionary = {}

# Directly requesting a list of all variants from the CIViC API in JSON format
variants = requests.get('https://civic.genome.wustl.edu/api/variants?count=10000').json()['records']

# Iterating through the entire list of variants and parsing out variants without coordinates, while keeping count of all variants
####################
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
####################
# Sorting all variants by stop coordinate in ascending order (least to greatest)
####################
raw_civic_variants_with_coordinates_list = sorted(raw_civic_variants_with_coordinates_list, key = lambda k: int(k['coordinates']['stop']))
for var in raw_civic_variants_with_coordinates_list:
	start = int(var['coordinates']['start'])
	stop = int(var['coordinates']['stop'])
	chr_key = var['coordinates']['chromosome']
	# Sorting nested and special cases of coordinates
	index_location = coordinate_sorter(start, stop, civic_variants_with_coordinates_dictionary[chr_key])
	# Inserting variants into civic_variants_with_coordinates_dictionary based on the index_location to maintain a sorted list
	civic_variants_with_coordinates_dictionary[chr_key].insert(index_location, var)
####################
##################################################
"""
for key,value in civic_variants_without_coordinates_dictionary.items():
	for index in range(0, len(value)):
		if value[index]['coordinates']['start'] == None:
			print(value[index]['entrez_name'] + '\t' + value[index]['name'])
"""
# This code read through everyline of an annotation file and maps those results back to CIViC
##################################################
found_mutations = ''
found_mutations_dict = {}

mutation_of_interest = ''
mutation_of_interest_dict = {}

#print(transvar_output.split('\n')[:5])

for line in transvar_output.split('\n'):
	#print(line)

	line_list = line.strip().split()

	#print(line_list)

	if line_list == []:
		break

	if line_list[0] == '#CHROM':

		for index in range(0, len(line_list)):
			if 'CHROM' in line_list[index]:
				chrom_index = index
			elif 'POS' in line_list[index]:
				start_index = index
			elif 'REF' in line_list[index]:
				ref_index = index
			elif 'ALT' in line_list[index]:
				alt_index = index
			elif 'transcript' in line_list[index]:
				rep_trans_index = index
			elif 'gene' in line_list[index]:
				gene_index = index+1
			elif 'coordinates' in line_list[index]:
				coordinates_index = index+1
			elif 'info' in line_list[index]:
				info_index = index+1

	if line_list[0].isdigit():

		var_ID = line_list[2]

		var_POS = line_list[start_index]

		var_REF = line_list[ref_index]

		var_ALT = line_list[alt_index]

		var_gene_name = line_list[gene_index]

		var_var_name = ''
		temp_name = line_list[coordinates_index].split('/')[-1]
		if 'p' in temp_name:
			var_var_name = line_list[coordinates_index].split('/')[-1].split('.')[1]
			var_var_short_name = var_var_name[:-1]
		elif 'fs' in temp_name:
			var_var_name = line_list[coordinates_index].split('/')[-1].split('.')[1]
			var_var_short_name = var_var_name.split('fs')[0][:-1]
		elif '.' == temp_name:
			var_var_name = None
			var_var_short_name = None
		
		var_chr = line_list[chrom_index]

		variant_start = ''
		variant_stop = ''
		if '_' in line_list[coordinates_index].split('/')[0].split('.')[1]:
			variant_start = line_list[info_index].split(';')[2].split('.')[1].split('_')[0]
			for index in range(0, len(line_list[info_index].split(';')[2].split('.')[1].split('_')[1])):
				if line_list[info_index].split(';')[2].split('.')[1].split('_')[1][index-1].isdigit() == True and line_list[info_index].split(';')[2].split('.')[1].split('_')[1][index].isalpha() == True:
					variant_stop = line_list[info_index].split(';')[2].split('.')[1].split('_')[1][:index]
		else:
			variant_start = line_list[start_index]
			for index in line_list[coordinates_index].split('/')[0].split('.')[1]:
				if index.isdigit() == True:
					variant_stop = variant_start

		temp_ref_base = len(line_list[ref_index])
		temp_alt_base = len(line_list[alt_index])
		if temp_ref_base > temp_alt_base:
			var_ref_base = line_list[ref_index][1:]
			var_var_base = None
		elif temp_ref_base < temp_alt_base:
			var_ref_base = None
			var_var_base = line_list[alt_index][2:]
		else:
			var_ref_base = line_list[ref_index]
			var_var_base = line_list[alt_index]

		var_rep_trans = line_list[rep_trans_index]

		# Implement genome build filtering once build 38 cooridnates are added to CIViC
		#genome_build = 37 OR 38


		#print(var_chr + '\t' + var_POS + '\t' + var_ID + '\t' + var_REF + '\t' + var_ALT + '\t' + var_gene_name + ':' + var_var_short_name + '\t' + var_chr + ':' + variant_start + '-' + variant_stop + '\t' + var_ref_base + '\t' + var_var_base + '\t' + var_rep_trans)

		civic_var_description = ''

		exact_match = ''
		passed_exact_match = False

		soft_match = ''
		passed_soft_match = False
		
		nested_match = ''
		nested_matchs = []
		passed_nested_match = False

		print_header_statement = False
		print_input_statement = False

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
					civic_var_description = current_variant_in_dict['description']

					found_mutations = var_chr + '\t' + var_POS + '\t' + var_ID + '\t' + var_REF + '\t' + var_ALT
					found_mutations_dict[found_mutations] = var_gene_name

					exact_match = current_variant_in_dict['entrez_name'] + ':' + current_variant_in_dict['name']

				## SOFT MATCHS ##
				elif (passed_exact_match == False and passed_soft_match == False) and ((var_var_short_name == civic_var_name) or (int(variant_start) in range(civic_start, civic_stop+1) and int(variant_stop) in range(civic_start, civic_stop+1)) or (civic_start in range(int(variant_start), int(variant_stop) + 1) and civic_stop in range(int(variant_start), int(variant_stop) + 1))):
				# elif (passed_exact_match == False and passed_soft_match == False and var_var_short_name == civic_var_name and var_ref_base != civic_ref_base and var_rep_trans == civic_rep_trans and var_gene_name == civic_gene_name) or 
				# 	(passed_exact_match == False and passed_soft_match == False and var_ref_base != civic_ref_base and int(variant_start) in range(civic_start, civic_stop+1) and int(variant_stop) in range(civic_start, civic_stop+1) and var_rep_trans == civic_rep_trans and var_gene_name == civic_gene_name) or 
				# 	(passed_exact_match == False and passed_soft_match == False and var_ref_base != civic_ref_base and civic_start in range(int(variant_start), int(variant_stop) + 1) and civic_stop in range(int(variant_start), int(variant_stop) + 1) and var_rep_trans == civic_rep_trans and var_gene_name == civic_gene_name) :
					passed_soft_match = True
					print_input_statement = True
					civic_var_description = current_variant_in_dict['description']

					found_mutations = var_chr + '\t' + var_POS + '\t' + var_ID + '\t' + var_REF + '\t' + var_ALT
					found_mutations_dict[found_mutations] = var_gene_name

					soft_match = current_variant_in_dict['entrez_name'] + ':' + current_variant_in_dict['name']

				## NESTED MATCHS ##
				elif (passed_exact_match == True or passed_soft_match == True) and var_var_name != civic_var_name and var_ref_base != civic_ref_base and int(variant_start) in range(civic_start, civic_stop+1) and int(variant_stop) in range( civic_start, civic_stop+1):
					passed_nested_match = True
					nested_match = current_variant_in_dict['entrez_name'] + ':' + current_variant_in_dict['name']
					nested_matchs.append(nested_match)

				## VARIANTS FOUND IN GENES IN CIViC BUT NOT MAPPED ##
				elif (passed_exact_match == False or passed_soft_match == False):
					mutation_of_interest = var_chr + '\t' + var_POS + '\t' + var_ID + '\t' + var_REF + '\t' + var_ALT
					mutation_of_interest_dict[mutation_of_interest] = var_gene_name
					#print(var_gene_name + ':.' + '\t' + var_chr + ':' + variant_start + '-' + variant_stop + '\t' + var_ref_base + '\t' + var_var_base + '\t' + var_rep_trans)

		# Printing mapping results to Stander Output
		if print_input_statement == True:
			if var_var_name == None:
				var_var_name = '.'
			if nested_matchs == []:
				nested_matchs.append('NONE')
			if civic_var_description == '':
				civic_var_description = 'No current description exist for this variant in CIViC'

				#print('#CHROM\tPOS\tID\tREF\tALT\tAnnotation\tMatch_level\tCIViC_variant_matched\tEncompassing_CIViC_variants\t')

			print('\nVariant input information:\tchr:' + var_chr + '\tStart:' + var_POS + '\tStop:' + variant_stop + '\tREF:' + var_REF + '\tALT:' + var_ALT)
			print('TransVar annotated variant as: ' + var_gene_name + ':' + var_var_name)
			print('\tMatches found to CIViC variants: ')
			
			if passed_exact_match == True:
				print '\tEXACT: ' + exact_match
			if passed_soft_match == True:
				print '\tSOFT: ' + soft_match
			if passed_nested_match == True:
				print '\t\tCIViC variants that might also be effected by variant:\n\t\t', u', '.join(nested_matchs)
			#print(var_chr + '\t' + var_POS + '\t' + var_ID + '\t' + var_REF + '\t' + var_ALT + '\t' + var_gene_name + ':' + var_var_name + '\t' + 'EXACT' + '\t' + exact_match + '\t' + str(nested_matchs).replace("'", ""))
				#print('\t' + civic_var_description.replace("\n", ""))

print('\nInputted variants occurring in genes found in CIViC but not mapped\n#CHROM\tPOS\tID\tREF\tALT\tCIViC_gene')
for key,value in mutation_of_interest_dict.items():
	if key not in found_mutations_dict:
		print(key + '\t' + value)
print('')




