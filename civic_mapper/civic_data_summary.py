
import argparse, json, requests, subprocess, transvar
import matplotlib.pyplot as plt
import numpy as np
requests.packages.urllib3.disable_warnings()


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

repeats = {}

publications = {}
unique_journals = {}
unique_drugs = {}
unique_disease = {}

evidence_items = requests.get('https://civic.genome.wustl.edu/api/evidence_items?count=10000000').json()['records']

raw_count = len(evidence_items)
for current_evidence in range(0, len(evidence_items)):

	# if evidence_items[current_evidence]["pubmed_id"] not in repeats and (evidence_items[current_evidence]["status"] == 'accepted' or evidence_items[current_evidence]["status"] == 'submitted'):
	# 	repeats[evidence_items[current_evidence]["pubmed_id"]] = [evidence_items[current_evidence]]
	# else:
	# 	repeats[evidence_items[current_evidence]["pubmed_id"]].append(evidence_items[current_evidence])


	if evidence_items[current_evidence]["pubmed_id"] not in publications:
		publications[evidence_items[current_evidence]["pubmed_id"]] = evidence_items[current_evidence]
		count += 1


# for key,value in repeats.items():
# 	status_list = []
# 	if len(value) > 1:
# 		print(key)
	# 	status_list.append(value[index]['status'])
	# if 'rejected' in status_list and 'accepted' in status_list:
	# 	print('\n\n' + value)



for key,value in publications.items():
	if value['status'] == 'rejected':
		#print('\n\n' + str(value))
		rejected += 1

	elif value['status'] == 'accepted':
		#print('\n\n' + str(value))
		accepted += 1

	elif value['status'] == 'submitted':
		#print(value['status'])
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
				#print(value['drugs'][index]['id'])
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
			#print(value['disease']['doid'])
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


#print('Raw evidence item count: ' + str(raw_count))
print('Unique pubmed IDs: ' + str(count))
print('\tAccepted: ' + str(accepted))
print('\tRejected: ' + str(rejected))
print('\tSubmitted: ' + str(submitted))
#print('\tTOTAL: ' + str(accepted+submitted))
fig = plt.figure()
plt.pie([accepted, rejected, submitted], labels = ['Accepted', 'Rejected', 'Submitted'], colors = ['red', 'gold', 'lightskyblue'], autopct='%1.1f%%', pctdistance = 0.65)
plt.axis('equal')
#plt.bar(left = [1,2,3] , height = [accepted, rejected, submitted], width = 0.8, align='center', tick_label = ['Accepted', 'Rejected', 'Submitted'])
plt.title('Total Publications (' + str(count) + ')')
fig.savefig('CIViC_publications.png')

# fig = plt.figure()
# plt.bar(left = [1,2,3], height = [accepted, rejected, submitted], align = 'center', color = ['red', 'gold', 'lightskyblue'], )
# plt.xticks([1,2,3], ('Accepted', 'Rejected', 'Submitted'))
# plt.title('Total Publications (' + str(count) + ')')
# fig.savefig('bar_CIViC_publications.png')


print('All stats from now on only take into account unique accepted or submitted journals with evidence items')
print('Unique Journals: ' + str(nonreject_art))
print('\tPrognostic: ' + str(prognostic))
print('\tDiagnostic: ' + str(diagnostic))
print('\tPredictive: ' + str(predictive))
#print('\tTOTAL: ' + str(prognostic+diagnostic+predictive))
fig = plt.figure()
plt.pie([prognostic, diagnostic, predictive], labels = ['Prognostic', 'Diagnostic', 'Predictive'], colors = ['red', 'gold', 'lightskyblue'], autopct='%1.1f%%', startangle = 90)
plt.axis('equal')
plt.title('Total Evidence Type (' + str(nonreject_art) + ')')
fig.savefig('CIViC_evidence_type.png')


print('\tClinical Significance:')
print('\t\tSensitivity: ' + str(sensitivity))
print('\t\tResistance or Non-Response: ' + str(resistance_non_response))
print('\t\tBetter Outcome: ' + str(better_outcome))
print('\t\tPoor Outcome: ' + str(poor_outcome))
print('\t\tPositive: ' + str(positive))
print('\t\tNegative: ' + str(negative))
print('\t\tAdverse Response: ' + str(adverse_response))
print('\t\tN/A: ' + str(n_a))
#print('\t\tTOTAL: ' + str(sensitivity+resistance_non_response+better_outcome+poor_outcome+positive+negative+adverse_response+n_a))
fig = plt.figure()
plt.pie([sensitivity, adverse_response, positive, resistance_non_response, negative, better_outcome, poor_outcome, n_a], labels = ['Sensitivity', 'Adverse Response', 'Positive', 'Resistance or Non-Response', 'Negative', 'Better Outcome', 'Poor Outcome', 'N/A'], autopct = '%1.1f%%', pctdistance=0.9, labeldistance=1.1, startangle = 40, colors = ['lightskyblue', 'white', 'yellowgreen', 'lightcyan', 'gold', 'green', 'lightcoral', 'red'])
plt.axis('equal')
plt.title('Total Clinical Significance (' + str(sensitivity+resistance_non_response+better_outcome+poor_outcome+positive+negative+adverse_response+n_a) + ')')
fig.savefig('CIViC_clinical_significance.png')

# fig = plt.figure()
# plt.bar(left = [1,2,3,4,5,6,7,8], height = [sensitivity, adverse_response, positive, resistance_non_response, negative, better_outcome, poor_outcome, n_a], align = 'center')
# plt.xticks([1,2,3,4,5,6,7,8], ('Sensitivity', 'Adverse Response', 'Positive', 'Resistance or Non-Response', 'Negative', 'Better Outcome', 'Poor Outcome', 'N/A'))
# plt.title('Total Clinical Significance (' + str(sensitivity+resistance_non_response+better_outcome+poor_outcome+positive+negative+adverse_response+n_a) + ')')
# fig.savefig('bar_CIViC_clinical_significance.png')


print('\tEvidence direction:')
print('\t\tSupports: ' + str(supports))
print('\t\tDoes Not Support: ' + str(does_not_support))
#print('\t\tTOTAL: ' + str(supports+does_not_support))
fig = plt.figure()
plt.pie([supports, does_not_support], labels = ['Supports', 'Does Not Support'], colors = ['gold', 'lightskyblue'], autopct='%1.1f%%', startangle = 180)
plt.axis('equal')
plt.title('Total Evidence Direction (' + str(supports+does_not_support) + ')')
fig.savefig('CIViC_evidence_direction.png')


print('\tTrust Ratings:')
print('\t\tTrust 1: ' + str(trust_rating_1))
print('\t\tTrust 2: ' + str(trust_rating_2))
print('\t\tTrust 3: ' + str(trust_rating_3))
print('\t\tTrust 4: ' + str(trust_rating_4))
print('\t\tTrust 5: ' + str(trust_rating_5))
#print('\t\tTOTAL: ' + str(trust_rating_1+trust_rating_2+trust_rating_3+trust_rating_4+trust_rating_5))
fig = plt.figure()
plt.pie([trust_rating_1, trust_rating_2, trust_rating_3, trust_rating_4, trust_rating_5], labels = ['1 Star', '2 Star', '3 Star', '4 Star', '5 Star'], autopct='%1.1f%%')
plt.axis('equal')
plt.title('Total Trust Ratings (' + str(trust_rating_1+trust_rating_2+trust_rating_3+trust_rating_4+trust_rating_5) + ')')
fig.savefig('CIViC_trust_ratings.png')


print('\tUnique Drugs: ' + str(drugs))
print('\tDiseases: ' + str(disease))
# fig = plt.figure()
# plt.bar(left = [1,2,3], height = [count, disease, drugs], align = 'center', color = ['red', 'gold', 'lightskyblue'], )
# plt.xticks([1,2,3], ('Publication', 'Diseases', 'Drugs'))
# plt.title('CIViC Statistics')
# fig.savefig('CIViC_statistics.png')


