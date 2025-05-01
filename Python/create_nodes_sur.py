# Import relevant modules
from py2neo import Graph


# Login to the database
graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")

query2 = """
        MATCH (file:File) WITH file

        WITH collect(file.url) AS fileURLs

        UNWIND fileURLs AS fileURL

        CALL apoc.periodic.iterate(

        '
        LOAD CSV WITH HEADERS FROM ($url) AS csv

        RETURN csv

        ',

        '

        WITH round(toFloat(csv.dendrogram_height),3,"HALF_UP") as height, csv.gene_name as gene, toInteger(csv.event_no) as event, csv.include as incl, round(toFloat(csv.median_logFC),3,"HALF_UP") as me,  toInteger(csv.leaf_status) as leaf, split(csv.patients, ",") as pts, toInteger(csv.no_of_uniq_patients) as no_pt, toInteger(csv.levels) as levels

        CREATE (:Expression {name: gene, event_no: event, dendrogram_height: height, median_exp: me, leaf_status: leaf, patients: pts, no_of_unique_patients: no_pt, include: incl, levels: levels});

        ',

        {batchSize:10000,parallel:true,params:{url:fileURL}}) YIELD batches, total

        RETURN batches, total;
        
        """
try:
    result2 = graph.run(query2)
    print(result2)
except:
    print('connection failed.............connect again!!!!!!!!!')
    graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")
    result2 = graph.run(query2)
    print(result2)

