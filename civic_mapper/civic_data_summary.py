"""This code queries the CIViC API for all evidence items,
filters out duplicated pubmed IDs, and returns counts of all
parameters listed along with pie chart summaries in pdf format"""

import argparse, json, requests, subprocess, transvar
import matplotlib.pyplot as plt

count = 0

rejected = 0
accepted = 0
submitted = 0

nonreject_art = 0

prognostic = 0
diagnostic = 0
predictive = 0

drugs = 0

sensitivity = 0
resistance_non_response = 0
better_outcome = 0
poor_outcome = 0
positive = 0
negative = 0
adverse_response = 0
n_a = 0

supports = 0
does_not_support = 0

disease = 0

trust_rating_1 = 0
trust_rating_2 = 0
trust_rating_3 = 0
trust_rating_4 = 0
trust_rating_5 = 0

publications = {}
unique_journals = {}
unique_drugs = {}
unique_disease = {}

evidence_items = requests.get('https://civic.genome.wustl.edu/api/evidence_items?count=10000000').json()['records']

for current_evidence in range(0, len(evidence_items)):
	if evidence_items[current_evidence]["pubmed_id"] not in publications:
		publications[evidence_items[current_evidence]["pubmed_id"]] = evidence_items[current_evidence]
		count += 1


for key,value in publications.items():
	if value['status'] == 'rejected':
		rejected += 1

	elif value['status'] == 'accepted':
		accepted += 1

	elif value['status'] == 'submitted':
		submitted += 1

	if value['status'] != 'rejected':
		if value['citation'].split(',')[-1].replace(".", "").upper() not in unique_journals:
			unique_journals[value['citation'].split(',')[-1].replace(".", "").upper()] = 1
			nonreject_art  += 1

		if value['evidence_type'] == 'Prognostic':
			prognostic += 1

		if value['evidence_type'] == 'Diagnostic':
			diagnostic += 1

		if value['evidence_type'] == 'Predictive':
			predictive += 1

		for index in range(0, len(value['drugs'])):
			if value['drugs'][index]['id'] not in unique_drugs:
				unique_drugs[value['drugs'][index]['id']] = 1
				drugs += 1

		if value['clinical_significance'] == 'Sensitivity':
			sensitivity += 1

		if value['clinical_significance'] == 'Resistance or Non-Response':
			resistance_non_response += 1

		if value['clinical_significance'] == 'Better Outcome':
			better_outcome += 1

		if value['clinical_significance'] == 'Poor Outcome':
			poor_outcome += 1

		if value['clinical_significance'] == 'Positive':
			positive += 1

		if value['clinical_significance'] == 'Negative':
			negative += 1

		if value['clinical_significance'] == 'Adverse Response':
			adverse_response += 1

		if value['clinical_significance'] == 'N/A':
			n_a += 1

		if value['evidence_direction'] == 'Supports':
			supports += 1

		if value['evidence_direction'] == 'Does Not Support':
			does_not_support += 1
		
		if value['disease']['doid'] not in unique_disease:
			unique_disease[value['disease']['doid']] = 1
			disease += 1

		if value['rating'] == 1:
			trust_rating_1 += 1

		if value['rating'] == 2:
			trust_rating_2 += 1
		
		if value['rating'] == 3:
			trust_rating_3 += 1

		if value['rating'] == 4:
			trust_rating_4 += 1

		if value['rating'] == 5:
			trust_rating_5 += 1

print('Unique Journals: ' + str(nonreject_art))
print('Diseases: ' + str(disease))
print('Unique Drugs: ' + str(drugs) + '\n')

print('Unique pubmed IDs: ' + str(count))
print('\tAccepted: ' + str(accepted))
print('\tRejected: ' + str(rejected))
print('\tSubmitted: ' + str(submitted) + '\n')
fig = plt.figure()
plt.pie([accepted, rejected, submitted], labels = ['Accepted', 'Rejected', 'Submitted'], colors = ['red', 'gold', 'lightskyblue'], autopct='%1.1f%%', pctdistance = 0.65)
plt.axis('equal')
plt.title('Total Publications (' + str(count) + ')')
fig.savefig('CIViC_publications.pdf', format = 'pdf')

print('All stats from now on only take into account unique accepted or submitted journals with evidence items')
print('\nEvidence Type Total: ' + str(prognostic+diagnostic+predictive))
print('\tPrognostic: ' + str(prognostic))
print('\tDiagnostic: ' + str(diagnostic))
print('\tPredictive: ' + str(predictive))
fig = plt.figure()
plt.pie([prognostic, diagnostic, predictive], labels = ['Prognostic', 'Diagnostic', 'Predictive'], colors = ['red', 'gold', 'lightskyblue'], autopct='%1.1f%%', startangle = 90)
plt.axis('equal')
plt.title('Total Evidence Type (' + str(prognostic+diagnostic+predictive) + ')')
fig.savefig('CIViC_evidence_type.pdf', format = 'pdf')

print('\nClinical Significance Total:' + str(sensitivity+resistance_non_response+better_outcome+poor_outcome+positive+negative+adverse_response+n_a))
print('\tSensitivity: ' + str(sensitivity))
print('\tResistance or Non-Response: ' + str(resistance_non_response))
print('\tBetter Outcome: ' + str(better_outcome))
print('\tPoor Outcome: ' + str(poor_outcome))
print('\tPositive: ' + str(positive))
print('\tNegative: ' + str(negative))
print('\tAdverse Response: ' + str(adverse_response))
print('\tN/A: ' + str(n_a))
fig = plt.figure()
plt.pie([sensitivity, adverse_response, positive, resistance_non_response, negative, better_outcome, poor_outcome, n_a], labels = ['Sensitivity', 'Adverse Response', 'Positive', 'Resistance or Non-Response', 'Negative', 'Better Outcome', 'Poor Outcome', 'N/A'], autopct = '%1.1f%%', pctdistance=0.9, labeldistance=1.1, startangle = 40, colors = ['lightskyblue', 'white', 'yellowgreen', 'lightcyan', 'gold', 'green', 'lightcoral', 'red'])
plt.axis('equal')
plt.title('Total Clinical Significance (' + str(sensitivity+resistance_non_response+better_outcome+poor_outcome+positive+negative+adverse_response+n_a) + ')')
fig.savefig('CIViC_clinical_significance.pdf', format = 'pdf')

print('\nEvidence direction Total: ' + str(supports+does_not_support))
print('\tSupports: ' + str(supports))
print('\tDoes Not Support: ' + str(does_not_support))
fig = plt.figure()
plt.pie([supports, does_not_support], labels = ['Supports', 'Does Not Support'], colors = ['gold', 'lightskyblue'], autopct='%1.1f%%', startangle = 180)
plt.axis('equal')
plt.title('Total Evidence Direction (' + str(supports+does_not_support) + ')')
fig.savefig('CIViC_evidence_direction.pdf', format = 'pdf')

print('\nTrust Ratings Total: ' + str(trust_rating_1+trust_rating_2+trust_rating_3+trust_rating_4+trust_rating_5))
print('\tTrust 1: ' + str(trust_rating_1))
print('\tTrust 2: ' + str(trust_rating_2))
print('\tTrust 3: ' + str(trust_rating_3))
print('\tTrust 4: ' + str(trust_rating_4))
print('\tTrust 5: ' + str(trust_rating_5))
fig = plt.figure()
plt.pie([trust_rating_1, trust_rating_2, trust_rating_3, trust_rating_4, trust_rating_5], labels = ['1 Star', '2 Star', '3 Star', '4 Star', '5 Star'], autopct='%1.1f%%')
plt.axis('equal')
plt.title('Total Trust Ratings (' + str(trust_rating_1+trust_rating_2+trust_rating_3+trust_rating_4+trust_rating_5) + ')')
fig.savefig('CIViC_trust_ratings.pdf', format = 'pdf')