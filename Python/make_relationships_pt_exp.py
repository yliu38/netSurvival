# Import relevant libraries
from py2neo import Graph

# -------------------------------------------------------------------------------------------------------------------#
# It is a good idea to store the filename into a variable.
# The variable can later become a function argument when the
# code is converted to a function body.

# connect to the server
graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")


query100 = """
        CALL apoc.periodic.iterate(

        '

        MATCH(n:Patient)

        MATCH(m:nodes_chose)
        
        WHERE n<>m AND any(x IN m.patients WHERE x IN n.name)

        RETURN id(n) as NId1, id(m) as NId2

        ',

        '

        MATCH (n:Patient), (m:nodes_chose)

        WHERE id(n) = NId1 AND id(m) = NId2

        CREATE (n)-[r:hasExpressionEvent]->(m)

        ',

        {batchSize:50000, parallel:true}) YIELD batches, total

        RETURN batches, total;
        """

try:
    result100 = graph.run(query100)
    print(result100)
except:
    print('connection failed.............connect again!!!!!!!!!')
    graph = Graph("neo4j://url", auth=("id", "password"), name = "dbname")
    result100 = graph.run(query100)
    print(result100)
