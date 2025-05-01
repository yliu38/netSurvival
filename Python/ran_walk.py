# Import relevant modules
from py2neo import Graph
import numpy as np
import pandas as pd
from time import time
import pickle
import multiprocessing as mp

# Kaplan-Meier curve
from sksurv.nonparametric import kaplan_meier_estimator
from lifelines.statistics import logrank_test

# Login to the database
graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")

def get_random_walks(k):

    # Query for random walk algorithm using the GDS library
    results = []
    # Length per walk
    for i in [3,6,9]:

        query = """
        CALL gds.randomWalk.stream(
            'ran_walk1',
                {
                walkLength: """  + str(i) + "," + \
                """walksPerNode: """ + str(k) + "," + \
                """randomSeed: 123,
                concurrency: 3
                }
        )
        YIELD nodeIds, path
        RETURN nodeIds, [node IN nodes(path) | node.name ] AS event_name
        """

        graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")
        # Run the random walk query
        result = graph.run(query)
        results.append(result.data())



    # Return the result
    flat_results = [x for sublist in results for x in sublist]
    return flat_results


def get_unique_patients(path):

    # Define a list to store all the patients in the path
    patient_list_path = []

    # For consecutive nodes in the random walk path, fetch the patients
    for first, second in zip(path, path[1:]):
        # print(first, second)

        # Query to fetch the patient ids for consecutive nodes
        query = """
        MATCH (p:Patient)--(n:nodes_chose) WHERE id(n) = """ + \
        str(first) + \
        """ WITH collect(p.name) as p1names, n MATCH (p:Patient)--(m:nodes_chose) WHERE id(m) = """ + \
        str(second) + \
        """ WITH p1names, n, collect(p.name) as p2names, m WITH n, m, apoc.coll.intersection(p1names, p2names) AS intersection RETURN intersection;"""

        # Run the query
        try:
            result = graph.run(query)
        except:
            print('connection failed.............connect again!!!!!!!!!')
            graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")
            result = graph.run(query)

        patient_list = result.data()[0]['intersection']
        # print(patient_list)

        # Concatinate the set of patients in the path
        patient_list_path = np.concatenate((patient_list_path, patient_list), axis = None)

    # Get the set of unique patients as a list
    unique_patients = list(set(list(patient_list_path)))

    # print(unique_patients)

    # Return unique patients
    return unique_patients

def get_all_patients():
    query = """
    MATCH (p:Patient) RETURN p.name
    """

    # Run the query
    try:
        result = graph.run(query)
    except:
        print('connection failed.............connect again!!!!!!!!!')
        graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")
        result = graph.run(query)

    # List to store all patients
    all_patients_list = []

    # Make a list to store all patients
    for x in result.data():
        all_patients_list.append(x['p.name'])

    return all_patients_list

def get_survival_data(unique_patients):

    # Make a list for getting the survival data
    survival_data = []
    censoring_data = []

    for x in unique_patients:
        graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")
        # Query survival time for each of those patients
        query = """MATCH (n:Patient) WHERE n.name = """ + \
        "'" + x + "'" + \
        """RETURN n.case_overall_survival_mo, n.censoring"""
        try:
            result = graph.run(query)
        except:
            print('connection failed.............connect again!!!!!!!!!')
            graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")
            result = graph.run(query)

        # Get the survival data from the query result
        tmp = result.data()[0]
        survival_data_patient = tmp['n.case_overall_survival_mo']
        survival_data = list(np.concatenate((survival_data, survival_data_patient), axis = None))

        # get censoring info
        censoring = tmp['n.censoring']
        censoring_data = list(np.concatenate((censoring_data, censoring), axis = None))

        columns = ['survival_data', 'censoring_data']
        tmp2 = pd.DataFrame(columns=columns, dtype=np.object_)
        tmp2['survival_data'] = survival_data
        tmp2['censoring_data'] = censoring_data

    # Get non-null survival data
    tmp3 = tmp2[tmp2.survival_data != 'NULL']
    tmp3['survival_data'] = tmp3['survival_data'].astype('float64')
    tmp3['censoring_data'] = tmp3['censoring_data'].astype('float64')
    return tmp3

def get_km_values(survival_data):

    #-----------------------------------------------------------------------------------#
    # KM survival analysis
    # Based on https://towardsdatascience.com/kaplan-meier-curve-explained-9c78a681faca
    #-----------------------------------------------------------------------------------#

    # Get the data corresponding to KM survival analysis
    # Get the duration and survival probabilities
    duration, survival_probability = kaplan_meier_estimator(survival_data['censoring_data']==1.0, survival_data['survival_data'])

    # Save them and return
    columns = ['duration', 'survival_probability']
    km = pd.DataFrame(columns=columns, dtype=np.float64)
    km['duration'] = duration
    km['survival_probability'] = survival_probability

    return km


def get_km_values_for_path(path):

    # Get unique patients in the path as a set
    unique_patients = get_unique_patients(path)

    # print(unique_patients)

    # Get the survival data for the unique patients as a list
    survival_data = get_survival_data(unique_patients)

    # Get Kaplan-Meier survival values
    km = get_km_values(survival_data)

    # print(km)
    # print(len(km))

    # Return km values
    return km, unique_patients, survival_data

def get_km_values_for_not_in_path(patients):

    # Get the survival data for the patients as a list
    survival_data = get_survival_data(patients)

    # Get Kaplan-Meier survival values
    km = get_km_values(survival_data)

    # Return km values
    return km, survival_data

def get_significant_paths(end_path):


    # Make a list for storing significant results
    sig_paths = []
    sig_kms = []
    sig_kms_not_in_path = []
    pvalues = []
    patients_in_each_path = []

    # Make an empty dictionary to store the results
    significant_results = {}

    # For each generated path evaluate survival for patients in the path
    # and not in the path
    for i in range(end_path):

        print("Working on path",i)

        # Get the path for the random walk i
        path = random_walks[i]['nodeIds']

        # Given a path this will get the KM values for that path and
        # the unique patients in that path
        print("...getting km values for the path")
        timeout = 180
        timeout_start = time()
        while True or time() < timeout_start + timeout:
            try:
                km, unique_patients, survival_data = get_km_values_for_path(path)
                break
            except:
                print('connection failed.............connect again!!!!!!!!!')
                graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")

        # Find the patients not in path
        patients_not_in_path = list(set(all_patients_list)-set(unique_patients))

        # Get the km values for patients not in path
        print("...getting km values for patients not in the path")
        timeout = 180
        timeout_start = time()
        while True or time() < timeout_start + timeout: 
            try:
                km_not, survival_not_data = get_km_values_for_not_in_path(patients_not_in_path)
                break
            except:
                print('connection failed.............connect again!!!!!!!!!')
                graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")


        # Use raw data and convert the values into float (does not consider censored info if using KM output)
        # km_duration_float = survival_data['survival_data'].apply(lambda x: float(x))
        # km_not_in_path_duration_float = survival_not_data['survival_data'].apply(lambda x: float(x))

        # Perform a log-rank test
        results = logrank_test(survival_data['survival_data'], survival_not_data['survival_data'], event_observed_A=survival_data['censoring_data'], event_observed_B=survival_not_data['censoring_data'])

        # Get the p-value
        pvalue = results.p_value
        pvalues.append(pvalue)

        # If the p-value is significant, store the correspinding paths, and associated data
        if pvalue < 0.05:

            print("...pvalue:", pvalue, "is **significant**, storing...")

            sig_paths.append(path)
            sig_kms.append(km)
            sig_kms_not_in_path.append(km_not)
            patients_in_each_path.append(unique_patients)
        else:
            print("...pvalue:", pvalue, "is not significant, moving on...")

    # Get the results as a dictionary
    significant_results = {'sig_paths':sig_paths,
                          'sig_kms':sig_kms,
                          'sig_kms_not_in_path':sig_kms_not_in_path,
                          'patients_in_each_path':patients_in_each_path}

    # Return all significant results
    return pvalues, significant_results


for xy in range(0,334):
    # First drop the graph if already exists
    query_drop = """RETURN gds.graph.exists('ran_walk1')"""
    try:
        result_drop = graph.run(query_drop)
        parsed_drop = result_drop.evaluate()
    except:
        print('connection failed.............connect again!!!!!!!!!')
        graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")
        result_drop = graph.run(query_drop)
        parsed_drop = result_drop.evaluate()

    if (parsed_drop == 1):
        try:
            result = graph.run("CALL gds.graph.drop('ran_walk1') YIELD graphName;")
            print("\nDropping query status for a graph named 'ran_walk1'\n")
        except:
            print('connection failed.............connect again!!!!!!!!!')
            graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")
            result = graph.run("CALL gds.graph.drop('ran_walk1') YIELD graphName;")
            print("\nDropping query status for a graph named 'ran_walk1'\n")


    # projection
    query1 = """
            // Cypher query
            MATCH (n:nodes_chose)
            WITH n.name AS geneName, rand() AS randomOrder, n
            ORDER BY geneName, randomOrder
            WITH geneName, COLLECT(n)[0] AS randomNode
            unwind(randomNode) as nme
            match (m:nodes_chose) WHERE id(m)=id(nme)
            set m:tmp1
            with m
            // Make a Cypher projection
            MATCH (source:tmp1)-[r:Top_Events]->(target:tmp1)
            WITH gds.graph.project('ran_walk1', source, target,{
            sourceNodeProperties: source {source: id(source) },
            targetNodeProperties: target {target: id(target) },
            relationshipProperties: r {num_common_patients: r.num_common_patients }}) AS g
            RETURN g.graphName AS graph, g.nodeCount AS nodes, g.relationshipCount AS rels
            """


    print("\nMade a subgraph and a Cypher projection named 'ran_walk1'\n")
    try:
        result1 = graph.run(query1)
        print(result1)
    except:
        print('connection failed.............connect again!!!!!!!!!')
        graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")
        result1 = graph.run(query1)
        print(result1)


    # random walks
    random_walks = get_random_walks(1)

    # Get the list of all patients
    all_patients_list = get_all_patients()

    # Define start and end paths for getting significant paths
    print(len(random_walks))

    # multiprocessing
    if __name__ == "__main__":
        start = time()
        p = mp.Pool(processes=8)
        result = p.map(get_significant_paths, [len(random_walks)])
        p.close()
        end = time()
        print("time taken = ", end-start)

    # Save the results

    with open('./pickle_rw/test_'+ str(xy) +'.pickle' , 'wb') as fl:
        pickle.dump(result, fl, pickle.HIGHEST_PROTOCOL)
        
       
    allpaths = [r['nodeIds'] for r in random_walks]
    with open('./pickle_rw/test_'+ str(xy) +'_allpath.pickle' , 'wb') as flp:
        pickle.dump(allpaths, flp, pickle.HIGHEST_PROTOCOL)
    
        
    query4 = """
            //delete tmp1
            match (n:tmp1)
            remove n:tmp1
            """
    print("\nRemove tmp label!\n")
    try:
        result4 = graph.run(query4)
        print(result4)
    except:
        print('connection failed.............connect again!!!!!!!!!')
        graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")
        result4 = graph.run(query4)
        print(result4)

