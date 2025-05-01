# -------------------------------------------------------------------------------------------------------------------#
# This python script will read the grouped cnv file and construct cql query for populating Neo4j database
# -------------------------------------------------------------------------------------------------------------------#
# Import relevant libraries
from py2neo import Graph

# -------------------------------------------------------------------------------------------------------------------#
# It is a good idea to store the filename into a variable.
# The variable can later become a function argument when the
# code is converted to a function body.
filename = '../pt_train.csv'

# connect to the server
graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")

# Using the newer with construct to close the file automatically.
with open(filename) as f:
    data = f.readlines()[1:]


for i in data:

    # Get the patient id
    pt_id = i.split(',')[6]
    
    # Get the sex
    case_sex = i.split(',')[1]
    
    # Get the case_vital_status
    case_vital_status = i.split(',')[3]
    
    # Get the case_overall_survival_mo
    case_overall_survival_mo = i.split(',')[4]

    # Get the censoring info
    censoring = i.split(',')[11]

    # Get the two year survival status
    status = i.split(',')[10]

    # Get the age
    age = i.split(',')[2]


 
    # add properties for existing nodes
    node = "MERGE (a:Patient" + "{name:" + str(pt_id)  + ", case_sex: " + str(case_sex) +", case_vital_status: " + str(case_vital_status) +", case_overall_survival_mo: " + "coalesce(" + str(case_overall_survival_mo) + ", 'NULL')"  + ", censoring: " + censoring + ", age: " + "coalesce(" + str(age) + ", 'NULL')"  + ", two_year_survival_status: " + status + "});"
           
    results = graph.run(node)
    print(results)

