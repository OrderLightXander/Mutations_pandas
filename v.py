import csv
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import click


iterator = 0
mutation_list = {}


def filterByCharacters(df):
    allowed_characters = set('ACTGN-')
    for p in df['Genome']:
        if not set(p).issubset(allowed_characters):
            df.drop((df[df['Genome'] == p ].index), inplace=True)
    return df

def detectMutation(dna, cdna):
    mutation_count = 1
    mutations = []
    for x, y in zip(dna, cdna):
        if x == y and x != '-' and y != '-':
            mutation_count += 1
        elif x != y and(x != 'N' or y !='N'): 
            #create mutation list for math plot
            if mutation_list.get(str(x)+str(mutation_count)+str(y)):
                mutation_list[str(x)+str(mutation_count)+str(y)] += 1
            else:
                mutation_list[str(x)+str(mutation_count)+str(y)] = 1
            #end instruction for plot
            mutations.append([x, mutation_count, y])
    return mutations, mutation_list

def selectGenomeByID(data, genome_id):
    genome = data.loc[data['ID'] == genome_id]
    return genome['Genome'].values[0]

def isDomainPresent(domain, data):
    if not re.findall('[^ACTGN]', domain):
        search_request = ''
        iterator = 0
        data['isDomainPresent'] = ''
        for character in domain:
            search_request += '[' + character + 'N]-*'
        for genome in data['Genome']:
            if re.findall(search_request, genome):
                data.at[iterator, 'isDomainPresent'] = True
                iterator+=1
            else :
                data.at[iterator, 'isDomainPresent'] = False
                iterator+=1
    else: 
        pass #Add desc. for non valid domain data
    
def buildVisualization(mutation_list, data, out_vis):
    sorted_mutations = dict(sorted(mutation_list.items(), key=lambda item: item[1], reverse=True))
    slised_mutations = dict(itertools.islice(sorted_mutations.items(), 10))
    for i in slised_mutations:
        slised_mutations[i] = slised_mutations[i]/len(data)
    f = plt.figure()
    plt.bar(*zip(*slised_mutations.items()))
    plt.show()
    f.savefig(out_vis, bbox_inches='tight')

@click.command()
@click.option('--input_csv', help='Number of greetings.') #'genomes.csv'
@click.option('--out_vis', help='Number of greetings.') #'plot.pdf'
@click.option('--out_csv', help='Number of greetings.') #'out.csv'
@click.option('--searched_domain', help='Number of greetings.') #'ACCCT'
@click.option('--reference_id', help='Number of greetings.') #'M_AURIS_581117'
def main(input_csv, out_vis, out_csv, searched_domain, reference_id):
    iterator = 0
    mutation_list = {}
    try:
        with open(input_csv, "r") as file:
            reader = csv.reader(file)
            df = pd.DataFrame(columns=['ID', 'Genome'])
            for row in reader:
                if row[0] != 'id':
                    df = df.append({'ID': row[0], 'Genome' : row[1]}, ignore_index=True)
    except IOError:
        print("File not acceseble")
    df = filterByCharacters(df).sort_values(by=['ID'])
    ds_sorted = df.reset_index(drop=True)
    ds_sorted['Mutations'] = ''
    ds_mutated = ds_sorted
    for z in ds_sorted['Genome']:
        if selectGenomeByID(ds_sorted, reference_id) != z:
            ds_mutated.at[iterator, 'Mutations'] = detectMutation(selectGenomeByID(ds_sorted, reference_id), z)[0]
            mutation_list.update(detectMutation(selectGenomeByID(ds_sorted, reference_id), z)[1])
            #print(detectMutation(selectGenomeByID(ds_sorted, reference_id), z)[1])
            iterator+=1
        else :
            ds_mutated.at[iterator, 'Mutations'] = 'Same DNA with referenced'
            iterator+=1
    ds_mutated.to_csv(out_csv) #'out.csv'
    isDomainPresent(searched_domain,ds_mutated)
    buildVisualization(mutation_list, ds_mutated, out_vis) #'plot.pdf'


if __name__ == '__main__':
    main()