#!/usr/bin/env python3

##########################################################
## For a given list of proteins the script resolves them
## (if possible) to the best matching STRING identifier
## and prints out the mapping on screen in the TSV format
##
## Requires requests module:
## type "python -m pip install requests" in command line
## (win) or terminal (mac/linux) to install the module
###########################################################

import requests ## python -m pip install requests
from time import sleep

####################################
## Format Genes for each Bacteria ##
####################################

## genes should be a list

# input = "./output/without_men/msea_significantRF/msea_results.csv" 
# output = "./output/without_men/string_significantRF"

input = "./output/without_men_saliva/msea_significantRF/msea_results.csv" 
output = "./output/without_men_saliva/string_significantRF"


gene_list = []
# bacteria_genes = {}
with open(input) as file:
    next(file)
    for line in file:

        ## use the q-value to only keep the significant gene interactions 
        if float(line.strip().split(",")[3]) < 0.049:
            # print( f"q-value: { float(line.strip().split(',')[3]) }" )
            gene = line.strip().split(",")[0]

            bacteria = line.strip().split("[\'")[1].split("\']")[0].replace("\'","").split(",")
            # if len(bacteria) == 3:
            if len(bacteria) == 1:
                gene_list.append(gene)

            # for bacterium in bacteria: 
            #     bacterium = bacterium.strip()
                
            #     if bacterium not in bacteria_genes:
            #         bacteria_genes[bacterium] = [gene]

            #     else:

            #         gene_list = bacteria_genes[bacterium] + [gene]
            #         bacteria_genes[bacterium] = gene_list

                

# print(bacteria_genes.keys())

# for bacteria, my_genes in bacteria_genes.items():    
print(gene_list)
############################
## Get STRING text output ##
############################
## get the string api url
string_api_url = "https://version-12-0.string-db.org/api"
output_format = "tsv"
method = "network"
# my_genes = list(my_genes)

## get the parameters
params = {

# "identifiers" : "\r".join(my_genes), # your protein list
"identifiers" : "%0d".join(gene_list), # the proteins
"species" : 9606, # NCBI/STRING taxon identifier: Human 
"add_white_nodes": 1,
# "limit" : 1, # only one (best) identifier per input protein
"network_flavor": "evidence", # show evidence links
"echo_query" : 1, # see your input identifiers in the output
"caller_identity" : "www.awesome_app.org", # your app name
"hide_disconnected_nodes": 1
}


## Construct URL
request_url = "/".join([string_api_url, output_format, method])

## Call STRING
results = requests.post(request_url, data=params)


## Read and parse the results
with open(f"{output}/allbacteria_network.tsv", "w") as file:
    for line in results.text.strip().split("\n"):
        file.write(f"{line}\n")
    


#############################
## Get STRING image output ##
#############################

string_api_url = "https://version-12-0.string-db.org/api"
output_format = "image"
method = "network"

##
## Construct URL
##


request_url = "/".join([string_api_url, output_format, method])

## For each gene call STRING
params = {

"identifiers" : "%0d".join(gene_list), # the proteins
"species" : 9606, # NCBI/STRING taxon identifier: Human 
"add_white_nodes": 0, # add 0 white nodes to my protein
"required_score": 0.9*1000, # same as minimum required interaction score on the website but on a different scale
"network_flavor": "evidence", # show evidence links
"caller_identity" : "www.awesome_app.org", # your app name
"hide_disconnected_nodes": 1 ## hide disconnected nodes
}


##
## Call STRING
##

response = requests.post(request_url, data=params)

##
## Save the network to file
##

file_name = f"{output}/allbacteria_network.png"
print("Saving interaction network to %s" % file_name)

with open(file_name, 'wb') as fh:
    fh.write(response.content)

sleep(1)







